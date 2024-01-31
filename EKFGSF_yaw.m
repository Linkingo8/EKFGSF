classdef EKFGSF_yaw < handle
    properties
        time = 0;
        doPrint = false;
        N_MODELS_EKFGSF = 5; % EKF实例数
        CONSTANTS_ONE_G = 9.80665; % m/s^2
        FLT_EPSILON = 1e-7;
        m_model_weights = 0;%1/N_MODELS_EKFGSF*ones(N_MODELS_EKFGSF,1);
        %% Parameters - these could be made tuneable
        m_gyro_noise = 0.1; 	% yaw rate noise used for covariance prediction (rad/sec)
        m_accel_noise = 2;		% horizontal accel noise used for covariance prediction (m/sec**2)
        m_tilt_gain = 0.2;		% gain from tilt error to gyro correction for complementary filter (1/sec)
        m_gyro_bias_gain = 0.04;	% gain applied to integral of gyro correction for complementary filter (1/sec)
        
        % Declarations used by the bank of N_MODELS_EKFGSF AHRS complementary filters
        m_true_airspeed = 0;	% true airspeed used for centripetal accel compensation (m/s)
        m_ahrs_ekf_gsf        
        
        m_ahrs_ekf_gsf_tilt_aligned = false;  % true the initial tilt alignment has been calculated
        m_ahrs_accel_fusion_gain = 0;      % gain from accel vector tilt error to rate gyro correction used by AHRS calculation
        m_ahrs_accel = zeros(3,1);     % low pass filtered body frame specific force vector used by AHRS calculation (m/s/s)
        m_ahrs_accel_norm = 0;             % length of _ahrs_accel specific force vector used by AHRS calculation (m/s/s)
        % Declarations used by a bank of N_MODELS_EKFGSF EKFs
        m_ekf_gsf
        m_vel_data_updated = 0;	% true when velocity data has been updated
        m_vel_NE = [0;0]; % NE velocity observations (m/s)
        m_vel_accuracy = 0;     % 1-sigma accuracy of velocity observations (m/s)
        m_ekf_gsf_vel_fuse_started = 0; % true when the EKF's have started fusing velocity data and the prediction and update processing is active
        
        % Declarations used by the Gaussian Sum Filter (GSF) that combines the individual EKF yaw estimates
        m_gsf_yaw = 0; 		% yaw estimate (rad)
        m_gsf_yaw_variance = 0; 	% variance of yaw estimate (rad^2)
    end
    methods
        function updateTime(obj, time)
            obj.time = time; 
        end
        function obj = EKFGSF_yaw(varargin)
            if length(varargin) == 1
                obj.time = varargin{1};
            end
            obj = initSL(obj);
            initialiseEKFGSF(obj, zeros(3,1), 0);
        end
        % Update Filter States - this should be called whenever new IMU data is available
        function obj = update(obj, ...
                imu_sample, ...
                run_EKF, ...  set to true when flying or movement is suitable for yaw estimation
                airspeed, ... true airspeed used for centripetal accel compensation - set to 0 when not required.
                imu_gyro_bias)  % estimated rate gyro bias (rad/sec)            
            % copy to class variables
            obj.m_true_airspeed = airspeed;
            
            % to reduce effect of vibration, filter using an LPF whose time constant is 1/10 of the AHRS tilt correction time constant
            filter_coef = min(10 * imu_sample.delta_vel_dt * obj.m_tilt_gain, 1);
            accel = imu_sample.delta_vel / max(imu_sample.delta_vel_dt, 0.001);
            obj.m_ahrs_accel = accel * (1 - filter_coef) + accel * filter_coef;
            
            % Initialise states first time
            if ~obj.m_ahrs_ekf_gsf_tilt_aligned
                % check for excessive acceleration to reduce likelihood of large initial roll/pitch errors
                % due to vehicle movement
                accel_norm_sq = norm_squared(obj, accel);
                accel_lpf_norm_sq = norm_squared(obj, obj.m_ahrs_accel);
                upper_accel_limit = obj.CONSTANTS_ONE_G * 1.1;
                lower_accel_limit = obj.CONSTANTS_ONE_G * 0.9;
                ok_to_align = ((accel_norm_sq > sq(obj, lower_accel_limit) && accel_norm_sq < sq(obj, upper_accel_limit))&& ...
                               (accel_lpf_norm_sq > sq(obj, lower_accel_limit) && accel_lpf_norm_sq < sq(obj, upper_accel_limit)));
                
                if (ok_to_align)
                    ahrsAlignTilt(obj, imu_sample.delta_vel);
                    obj.m_ahrs_ekf_gsf_tilt_aligned = true;
                end
                
                return;
            end
            
            % calculate common values used by the AHRS complementary filter models
            obj.m_ahrs_accel_norm = norm(obj.m_ahrs_accel);
            
            % AHRS prediction cycle for each model - this always runs
            obj.m_ahrs_accel_fusion_gain = obj.ahrsCalcAccelGain();
                    
            for model_index = 1:obj.N_MODELS_EKFGSF
                predictEKF(obj, model_index, imu_sample.delta_ang, imu_sample.delta_ang_dt, imu_sample.delta_vel, imu_sample.delta_vel_dt);
            end
            
            % The 3-state EKF models only run when flying to avoid corrupted estimates due to operator handling and GPS interference
            if (run_EKF && obj.m_vel_data_updated)
                if ~obj.m_ekf_gsf_vel_fuse_started
                    initialiseEKFGSF(obj, obj.m_vel_NE, obj.m_vel_accuracy);
                    
                    % Initialise to gyro bias estimate from main filter because there could be a large
                    % uncorrected rate gyro bias error about the gravity vector
                    ahrsAlignYaw(obj, imu_gyro_bias);
                    
                    if(sqrt(sq(obj.m_vel_NE(1))+sq(obj.m_vel_NE(2)))>abs(obj.m_vel_accuracy))
                        obj.m_ekf_gsf_vel_fuse_started = true;
                    end
                    obj.m_vel_data_updated = false;
                else
                    bad_update = false;
                    for model_index = 1:obj.N_MODELS_EKFGSF
                        % subsequent measurements are fused as direct state observations                        
                        if ~updateEKF(obj, model_index, obj.m_vel_NE, obj.m_vel_accuracy)
                            bad_update = true;
                        end
                    end
                    
                    if ~bad_update
                        total_weight = 0;
                        % calculate weighting for each model assuming a normal distribution
                        min_weight = 1e-5;
                        n_weight_clips = 0;
                        for model_index = 1:obj.N_MODELS_EKFGSF
                            obj.m_model_weights(model_index) = gaussianDensity(obj, model_index) * obj.m_model_weights(model_index);
                            if obj.m_model_weights(model_index) < min_weight
                                n_weight_clips = n_weight_clips + 1;
                                obj.m_model_weights(model_index) = min_weight;
                            end
                            
                            total_weight = total_weight + obj.m_model_weights(model_index);
                        end
                        
                        % normalise the weighting function             
                        if n_weight_clips < obj.N_MODELS_EKFGSF
                            obj.m_model_weights = obj.m_model_weights/total_weight;
                        else
                            % all weights have collapsed due to excessive innovation variances so reset filters
                            fprintf('reset filters\n');
                            initialiseEKFGSF(obj, obj.m_vel_NE, obj.m_vel_accuracy);
                        end
                    end

                    % Calculate a composite yaw vector as a weighted average of the states for each model.
                    % To avoid issues with angle wrapping, the yaw state is converted to a vector with length
                    % equal to the weighting value before it is summed.
                    yaw_vector = [0;0];
                    
                    for model_index = 1:obj.N_MODELS_EKFGSF
                        yaw_vector(1) = yaw_vector(1) + obj.m_model_weights(model_index) * cos(obj.m_ekf_gsf(model_index).X(3));
                        yaw_vector(2) = yaw_vector(2) + obj.m_model_weights(model_index) * sin(obj.m_ekf_gsf(model_index).X(3));
                    end
                    
                    obj.m_gsf_yaw = atan2(yaw_vector(2), yaw_vector(1));
                    
                    % calculate a composite variance for the yaw state from a weighted average of the variance for each model
                    % models with larger innovations are weighted less
                    obj.m_gsf_yaw_variance = 0;
                    
                    for model_index = 1:obj.N_MODELS_EKFGSF
                        yaw_delta = wrap_pi(obj, obj.m_ekf_gsf(model_index).X(3) - obj.m_gsf_yaw);
                        obj.m_gsf_yaw_variance = obj.m_gsf_yaw_variance + obj.m_model_weights(model_index) * (obj.m_ekf_gsf(model_index).P(3, 3) + yaw_delta * yaw_delta);
                    end
                    
                    % prevent the same velocity data being used more than once
                    obj.m_vel_data_updated = false;
                end
                
            elseif obj.m_ekf_gsf_vel_fuse_started && ~run_EKF
                % wait to fly again
                obj.m_ekf_gsf_vel_fuse_started = false;
            end
           
        end
        function obj = setVelocity(obj, velocity, ... NE velocity measurement (m/s)
                accuracy)	% 1-sigma accuracy of velocity measurement (m/s)
            obj.m_vel_NE = velocity;
            obj.m_vel_accuracy = accuracy;
            obj.m_vel_data_updated = true;
        end

        function obj = setGyroBias(obj,imu_gyro_bias)
              for model_index = 1:obj.N_MODELS_EKFGSF
                    obj.m_ahrs_ekf_gsf(model_index).gyro_bias = imu_gyro_bias;
              end
        end
        % get solution data for logging
        function out = getLogData(obj, yaw_composite,...
                yaw_composite_variance,...
                yaw,...      维数 N_MODELS_EKFGSF
                innov_VN,... 维数 N_MODELS_EKFGSF
                innov_VE,... 维数 N_MODELS_EKFGSF
                weight)

        end
        function out = isActive(obj)
            out = obj.m_ekf_gsf_vel_fuse_started;
        end
        function  out = getYaw(obj)
            out = obj.m_gsf_yaw;
        end
        function out = getYawVar(obj)
            out = obj.m_gsf_yaw_variance;
        end
        
        
        % calculate the gain from gravity vector misalingment to tilt correction to be used by all AHRS filters
        function out = ahrsCalcAccelGain(obj)
            % Calculate the acceleration fusion gain using a continuous function that is unity at 1g and zero
            % at the min and max g value. Allow for more acceleration when flying as a fixed wing vehicle using centripetal
            % acceleration correction as higher and more sustained g will be experienced.
            % Use a quadratic instead of linear function to prevent vibration around 1g reducing the tilt correction effectiveness.
            % see https://www.desmos.com/calculator/dbqbxvnwfg
            
            attenuation = 2;
            centripetal_accel_compensation_enabled = obj.m_true_airspeed > obj.FLT_EPSILON;
            
            if centripetal_accel_compensation_enabled && obj.m_ahrs_accel_norm > obj.CONSTANTS_ONE_G
                attenuation = 1;
            end
            
            delta_accel_g = (obj.m_ahrs_accel_norm - obj.CONSTANTS_ONE_G) / obj.CONSTANTS_ONE_G;
            out = obj.m_tilt_gain * obj.sq(1 - min(attenuation * abs(delta_accel_g), 1));
        end
        %update specified AHRS rotation matrix using IMU and optionally true airspeed data
        function ahrsPredict(obj, model_index, delta_ang, delta_ang_dt)
            % generate attitude solution using simple complementary filter for the selected model
            ang_rate = delta_ang / max(delta_ang_dt, 0.001) - obj.m_ahrs_ekf_gsf(model_index).gyro_bias;
            
            R_to_body = obj.m_ahrs_ekf_gsf(model_index).R';
            gravity_direction_bf = R_to_body(:,3);
            
            % Perform angular rate correction using accel data and reduce correction as accel magnitude moves away from 1 g (reduces drift when vehicle picked up and moved).
            % During fixed wing flight, compensate for centripetal acceleration assuming coordinated turns and X axis forward
            tilt_correction = zeros(3,1);
            
            if obj.m_ahrs_accel_fusion_gain > 0.0
                
                accel = obj.m_ahrs_accel;
                
                if obj.m_true_airspeed > obj.FLT_EPSILON
                    % Calculate body frame centripetal acceleration with assumption X axis is aligned with the airspeed vector
                    % Use cross product of body rate and body frame airspeed vector
                    centripetal_accel_bf = [0; obj.m_true_airspeed * ang_rate(3); - obj.m_true_airspeed * ang_rate(2)];
                    
                    % correct measured accel for centripetal acceleration
                    accel = accel - centripetal_accel_bf;
                end
                
                tilt_correction = cross(gravity_direction_bf , accel) * obj.m_ahrs_accel_fusion_gain / obj.m_ahrs_accel_norm;
            end
            
            % Gyro bias estimation
            gyro_bias_limit = 0.05;
            spin_rate = norm(ang_rate);
            
            if spin_rate < 10*pi/180                                                                                        %%误差积分项
                obj.m_ahrs_ekf_gsf(model_index).gyro_bias = obj.m_ahrs_ekf_gsf(model_index).gyro_bias - tilt_correction * (obj.m_gyro_bias_gain * delta_ang_dt);
                obj.m_ahrs_ekf_gsf(model_index).gyro_bias = min(max(obj.m_ahrs_ekf_gsf(model_index).gyro_bias,-gyro_bias_limit), gyro_bias_limit);
            else
%                 fprintf("spin rate : %f to large \n",spin_rate);
            end
            
            % delta angle from previous to current frame
            delta_angle_corrected = delta_ang + (tilt_correction - obj.m_ahrs_ekf_gsf(model_index).gyro_bias) * delta_ang_dt;
            % Apply delta angle to rotation matrix
            obj.m_ahrs_ekf_gsf(model_index).R = obj.ahrsPredictRotMat(obj.m_ahrs_ekf_gsf(model_index).R, delta_angle_corrected);
        end
        % align all AHRS roll and pitch orientations using IMU delta velocity vector
        function ahrsAlignTilt(obj, delta_vel)
            % Rotation matrix is constructed directly from acceleration measurement and will be the same for
            % all models so only need to calculate it once. Assumptions are:
            % 1) Yaw angle is zero - yaw is aligned later for each model when velocity fusion commences.
            % 2) The vehicle is not accelerating so all of the measured acceleration is due to gravity.
            if obj.doPrint
                fprintf('ahrsAlignTilt %.3f\n',obj.time);
            end
            % Calculate earth frame Down axis unit vector rotated into body frame
            down_in_bf = -obj.normalized(delta_vel);
            
            % Calculate earth frame North axis unit vector rotated into body frame, orthogonal to 'down_in_bf'
            i_vec_bf = [1;0;0];
            north_in_bf = i_vec_bf - down_in_bf * dot(i_vec_bf,down_in_bf);
            north_in_bf = obj.normalized(north_in_bf);
            
            % Calculate earth frame East axis unit vector rotated into body frame, orthogonal to 'down_in_bf' and 'north_in_bf'
            east_in_bf = cross(down_in_bf , north_in_bf);
            
            % Each column in a rotation matrix from earth frame to body frame represents the projection of the
            % corresponding earth frame unit vector rotated into the body frame, eg 'north_in_bf' would be the first column.
            % We need the rotation matrix from body frame to earth frame so the earth frame unit vectors rotated into body
            % frame are copied into corresponding rows instead.
            R = [north_in_bf';east_in_bf';down_in_bf'];
            
            for model_index = 1:obj.N_MODELS_EKFGSF
                obj.m_ahrs_ekf_gsf(model_index).R = R;
            end
        end
        %align all AHRS yaw orientations to initial values
        function ahrsAlignYaw(obj, imu_gyro_bias)
            % Align yaw angle for each model
            for model_index = 1:obj.N_MODELS_EKFGSF
                R = obj.m_ahrs_ekf_gsf(model_index).R;
                yaw = obj.wrap_pi(obj.m_ekf_gsf(model_index).X(3));
                obj.m_ahrs_ekf_gsf(model_index).R = obj.updateYawInRotMat(yaw, R);
                
%                 obj.m_ahrs_ekf_gsf(model_index).aligned = true;
                
%                 obj.m_ahrs_ekf_gsf(model_index).gyro_bias = imu_gyro_bias;
            end
        end
        % Efficient propagation of a delta angle in body frame applied to the body to earth frame rotation matrix
        function ret = ahrsPredictRotMat(obj, R, g)
            ret = R;
            ret(1, 1) = ret(1, 1) + R(1, 2) * g(3) - R(1, 3) * g(2);
            ret(1, 2) = ret(1, 2) + R(1, 3) * g(1) - R(1, 1) * g(3);
            ret(1, 3) = ret(1, 3) + R(1, 1) * g(2) - R(1, 2) * g(1);
            ret(2, 1) = ret(2, 1) + R(2, 2) * g(3) - R(2, 3) * g(2);
            ret(2, 2) = ret(2, 2) + R(2, 3) * g(1) - R(2, 1) * g(3);
            ret(2, 3) = ret(2, 3) + R(2, 1) * g(2) - R(2, 2) * g(1);
            ret(3, 1) = ret(3, 1) + R(3, 2) * g(3) - R(3, 3) * g(2);
            ret(3, 2) = ret(3, 2) + R(3, 3) * g(1) - R(3, 1) * g(3);
            ret(3, 3) = ret(3, 3) + R(3, 1) * g(2) - R(3, 2) * g(1);
            
            % Renormalise rows
            for r = 1:3
                rowLengthSq = obj.norm_squared(ret(r,:));
                
                if (rowLengthSq > obj.FLT_EPSILON)
                    % Use linear approximation for inverse sqrt taking advantage of the row length being close to 1.0
                    rowLengthInv = 1.5 - 0.5 * rowLengthSq;
                    ret(r,:) = ret(r,:) * rowLengthInv;
                end
            end
        end
        
        function obj = initSL(obj)
            obj.m_ahrs_ekf_gsf(obj.N_MODELS_EKFGSF).R = eye(3);             % matrix that rotates a vector from body to earth frame
            obj.m_ahrs_ekf_gsf(obj.N_MODELS_EKFGSF).gyro_bias = zeros(3,1); % gyro bias learned and used by the quaternion calculation
            obj.m_ahrs_ekf_gsf(obj.N_MODELS_EKFGSF).aligned = false;        % true when AHRS has been aligned
            for i = 1:obj.N_MODELS_EKFGSF
                obj.m_ahrs_ekf_gsf(i).R = eye(3);             % matrix that rotates a vector from body to earth frame
                obj.m_ahrs_ekf_gsf(i).gyro_bias = zeros(3,1); % gyro bias learned and used by the quaternion calculation
                obj.m_ahrs_ekf_gsf(i).aligned = false;        % true when AHRS has been aligned
            end
            
        end
        function initialiseEKFGSF(obj, vel_NE, vel_accuracy)
            if obj.doPrint
                fprintf('initialiseEKFGSF %.3f\n',obj.time);
            end            
            obj.m_gsf_yaw = 0;
            obj.m_ekf_gsf_vel_fuse_started = false;
            obj.m_gsf_yaw_variance = obj.sq(pi / 2);
            obj.m_model_weights = 1/obj.N_MODELS_EKFGSF*ones(obj.N_MODELS_EKFGSF,1);  % All filter models start with the same weight
            
            obj.m_ekf_gsf(obj.N_MODELS_EKFGSF).X = zeros(3,1);         % Vel North (m/s),  Vel East (m/s), yaw (rad)s
            obj.m_ekf_gsf(obj.N_MODELS_EKFGSF).P = zeros(3,3);         % covariance matrix
            obj.m_ekf_gsf(obj.N_MODELS_EKFGSF).S_inverse = zeros(2,2); % inverse of the innovation covariance matrix
            obj.m_ekf_gsf(obj.N_MODELS_EKFGSF).S_det_inverse = 0;      % inverse of the innovation covariance matrix determinant
            obj.m_ekf_gsf(obj.N_MODELS_EKFGSF).innov = zeros(2,1);     % Velocity N,E innovation (m/s)5
            for i = 1:obj.N_MODELS_EKFGSF
                obj.m_ekf_gsf(i).X = zeros(3,1);         % Vel North (m/s),  Vel East (m/s), yaw (rad)s
                obj.m_ekf_gsf(i).P = zeros(3,3);         % covariance matrix
                obj.m_ekf_gsf(i).S_inverse = zeros(2,2); % inverse of the innovation covariance matrix
                obj.m_ekf_gsf(i).S_det_inverse = 0;      % inverse of the innovation covariance matrix determinant
                obj.m_ekf_gsf(i).innov = zeros(2,1);     % Velocity N,E innovation (m/s)5
            end
            yaw_increment = 2 * pi / obj.N_MODELS_EKFGSF;
            for model_index = 1:obj.N_MODELS_EKFGSF
                % evenly space initial yaw estimates in the region between +-Pi
                obj.m_ekf_gsf(model_index).X(3) = -pi + 0.5 * (yaw_increment) + (model_index-1) * yaw_increment;
                % take velocity states and corresponding variance from last measurement
                obj.m_ekf_gsf(model_index).X(1) = vel_NE(1);
                obj.m_ekf_gsf(model_index).X(2) = vel_NE(2);
                obj.m_ekf_gsf(model_index).P(1, 1) = obj.sq(max(vel_accuracy, 0.01));
                obj.m_ekf_gsf(model_index).P(2, 2) = obj.m_ekf_gsf(model_index).P(1, 1);
                % use half yaw interval for yaw uncertainty
                obj.m_ekf_gsf(model_index).P(3, 3) = obj.sq(0.5 * yaw_increment);
            end
        end
        % predict state and covariance for the specified EKF using inertial data
        % delta_ang 3x1
        % delta_ang_dt 1x1
        function predictEKF(obj, model_index, delta_ang, delta_ang_dt, delta_vel, delta_vel_dt)
            % generate an attitude reference using IMU data
            obj.ahrsPredict(model_index, delta_ang, delta_ang_dt);
            
            % we don't start running the EKF part of the algorithm until there are regular velocity observations
            if (~obj.m_ekf_gsf_vel_fuse_started)
                return;
            end
            
            % Calculate the yaw state using a projection onto the horizontal that avoids gimbal lock
            obj.m_ekf_gsf(model_index).X(3) = getEulerYaw(obj, obj.m_ahrs_ekf_gsf(model_index).R);
            
            % calculate delta velocity in a horizontal front-right frame
            del_vel_NED = obj.m_ahrs_ekf_gsf(model_index).R * delta_vel;
            cos_yaw = cos(obj.m_ekf_gsf(model_index).X(3));
            sin_yaw = sin(obj.m_ekf_gsf(model_index).X(3));
            dvx =   del_vel_NED(1) * cos_yaw + del_vel_NED(2) * sin_yaw;
            dvy = - del_vel_NED(1) * sin_yaw + del_vel_NED(2) * cos_yaw;
            
            del_ang_NED = obj.m_ahrs_ekf_gsf(model_index).R * delta_ang;
	        daz = del_ang_NED(3);

            % Use fixed values for delta velocity and delta angle process noise variances
            d_vel_var = obj.sq(2*obj.m_accel_noise * delta_vel_dt);
            d_ang_var = obj.sq(obj.m_gyro_noise * delta_ang_dt);
            
            obj.m_ekf_gsf(model_index).P = YawEstPredictCovariance(obj, obj.m_ekf_gsf(model_index).X, obj.m_ekf_gsf(model_index).P, [dvx; dvy],d_vel_var, daz,d_ang_var);
            
            % covariance matrix is symmetrical, so copy upper half to lower half
            obj.m_ekf_gsf(model_index).P(2, 1) = obj.m_ekf_gsf(model_index).P(1, 2);
            obj.m_ekf_gsf(model_index).P(3, 1) = obj.m_ekf_gsf(model_index).P(1, 3);
            obj.m_ekf_gsf(model_index).P(3, 2) = obj.m_ekf_gsf(model_index).P(2, 3);
            
            % constrain variances
            min_var = 1e-6;
            
            for index = 1:3
                obj.m_ekf_gsf(model_index).P(index, index) = max(obj.m_ekf_gsf(model_index).P(index, index), min_var);
            end
            
            % sum delta velocities in earth frame:
            obj.m_ekf_gsf(model_index).X(1) = obj.m_ekf_gsf(model_index).X(1) + del_vel_NED(1);
            obj.m_ekf_gsf(model_index).X(2) = obj.m_ekf_gsf(model_index).X(2) + del_vel_NED(2);
        end
        % update state and covariance for the specified EKF using a NE velocity measurement
        % return false if update failed
        function out = updateEKF(obj, model_index, vel_NE, vel_accuracy)

            % set observation variance from accuracy estimate supplied by GPS and apply a sanity check minimum
            vel_obs_var = obj.sq(max(vel_accuracy, 0.01));
            
            % calculate velocity observation innovations
            obj.m_ekf_gsf(model_index).innov(1) = obj.m_ekf_gsf(model_index).X(1) - vel_NE(1);
            obj.m_ekf_gsf(model_index).innov(2) = obj.m_ekf_gsf(model_index).X(2) - vel_NE(2);
            
            K = zeros(3,2);
            P_new = zeros(3,3);%	matrix::SquareMatrix<float, 3> P_new;
            
            [obj, S_inv, S_det_inv, K, P_new] = YawEstComputeMeasurementUpdate(obj, obj.m_ekf_gsf(model_index).P, vel_obs_var, obj.FLT_EPSILON);
            obj.m_ekf_gsf(model_index).S_inverse = S_inv;
            obj.m_ekf_gsf(model_index).S_det_inverse = S_det_inv;
            
            obj.m_ekf_gsf(model_index).P = P_new;
            
            % copy upper to lower diagonal
            obj.m_ekf_gsf(model_index).P(2, 1) = obj.m_ekf_gsf(model_index).P(1, 2);
            obj.m_ekf_gsf(model_index).P(3, 1) = obj.m_ekf_gsf(model_index).P(1, 3);
            obj.m_ekf_gsf(model_index).P(3, 2) = obj.m_ekf_gsf(model_index).P(2, 3);
            
            % constrain variances
            min_var = 1e-6;
            
            for index = 1:3
                obj.m_ekf_gsf(model_index).P(index, index) = max(obj.m_ekf_gsf(model_index).P(index, index), min_var);
            end
            
            % test ratio = transpose(innovation) * inverse(innovation variance) * innovation = [1x2] * [2,2] * [2,1] = [1,1]
            test_ratio = obj.m_ekf_gsf(model_index).innov' * (obj.m_ekf_gsf(model_index).S_inverse * obj.m_ekf_gsf(model_index).innov);
            
            % Perform a chi-square innovation consistency test and calculate a compression scale factor
            % that limits the magnitude of innovations to 5-sigma
            % If the test ratio is greater than 25 (5 Sigma) then reduce the length of the innovation vector to clip it at 5-Sigma
            % This protects from large measurement spikes
            if test_ratio > 25
                innov_comp_scale_factor = sqrt(25 / test_ratio);
            else
                innov_comp_scale_factor = 1;
            end
            % Correct the state vector and capture the change in yaw angle
            oldYaw = obj.m_ekf_gsf(model_index).X(3);
            
            obj.m_ekf_gsf(model_index).X = obj.m_ekf_gsf(model_index).X - (K * obj.m_ekf_gsf(model_index).innov) * innov_comp_scale_factor;
            
            yawDelta = obj.m_ekf_gsf(model_index).X(3) - oldYaw;
            
            % apply the change in yaw angle to the AHRS
            % take advantage of sparseness in the yaw rotation matrix
            cosYaw = cos(yawDelta);
            sinYaw = sin(yawDelta);
            R_prev00 = obj.m_ahrs_ekf_gsf(model_index).R(1, 1);
            R_prev01 = obj.m_ahrs_ekf_gsf(model_index).R(1, 2);
            R_prev02 = obj.m_ahrs_ekf_gsf(model_index).R(1, 3);
            
            obj.m_ahrs_ekf_gsf(model_index).R(1, 1) = R_prev00 * cosYaw - obj.m_ahrs_ekf_gsf(model_index).R(2, 1) * sinYaw;
            obj.m_ahrs_ekf_gsf(model_index).R(1, 2) = R_prev01 * cosYaw - obj.m_ahrs_ekf_gsf(model_index).R(2, 2) * sinYaw;
            obj.m_ahrs_ekf_gsf(model_index).R(1, 3) = R_prev02 * cosYaw - obj.m_ahrs_ekf_gsf(model_index).R(2, 3) * sinYaw;
            obj.m_ahrs_ekf_gsf(model_index).R(2, 1) = R_prev00 * sinYaw + obj.m_ahrs_ekf_gsf(model_index).R(2, 1) * cosYaw;
            obj.m_ahrs_ekf_gsf(model_index).R(2, 2) = R_prev01 * sinYaw + obj.m_ahrs_ekf_gsf(model_index).R(2, 2) * cosYaw;
            obj.m_ahrs_ekf_gsf(model_index).R(2, 3) = R_prev02 * sinYaw + obj.m_ahrs_ekf_gsf(model_index).R(2, 3) * cosYaw;
            
            out = true;
        end
        function out = sq(obj,x)
            out = x^2;
        end
        
        function out = norm_squared(obj,x)
            out = dot(x,x);
        end
        
        function out = wrap_pi(obj, x)
            n_angle = rem(x, 2 * pi);
        
            if n_angle > pi
                out = n_angle - 2*pi;
            elseif n_angle < -pi
                out = n_angle + 2*pi;
            else 
                out = n_angle;
            end
        end
        
        function out = normalized(obj,x)
            out = x/norm(x);
        end
        % return the probability of the state estimate for the specified EKF assuming a gaussian error distribution
        function out = gaussianDensity(obj, model_index)
            normDist = dot(obj.m_ekf_gsf(model_index).innov, obj.m_ekf_gsf(model_index).S_inverse * obj.m_ekf_gsf(model_index).innov);
            out = 1/(2 * pi) * sqrt(obj.m_ekf_gsf(model_index).S_det_inverse) * exp(-0.5 * normDist);           
        end
        % Checks which euler rotation sequence to use and update yaw in rotation matrix
        function out = updateYawInRotMat(obj, yaw, rot_in)
            if obj.shouldUse321RotationSequence(rot_in)
                out = obj.updateEuler321YawInRotMat(yaw, rot_in);
%                 fprintf('321Rotation\n');
            else
                out = obj.updateEuler312YawInRotMat(yaw, rot_in);
%                 fprintf('312Rotation\n');
            end
        end
        %
        function out = getEulerYaw(obj, R)
            if obj.shouldUse321RotationSequence(R)
                out = obj.getEuler321Yaw(R);
%                 fprintf('321Rotation\n');
            else
                out = obj.getEuler312Yaw(R);
%                 fprintf('312Rotation\n');
            end
        end
        function out = getEuler321Yaw(obj, R)
            out = atan2(R(2, 1), R(1, 1));
        end
        function out = getEuler312Yaw(obj, R)
            out = atan2(-R(1, 2), R(2, 2));
        end
        %
        function out = shouldUse321RotationSequence(obj,R)
            out = false;
            if abs(R(3,1)) < abs(R(3,2))  % [index]
                out = true;
            end
        end
        %
        function out = updateEuler312YawInRotMat(obj, yaw, rot_in)
            rotVec312 = [yaw; asin(rot_in(3, 2)); atan2(-rot_in(3, 1), rot_in(3, 3))]; %yaw roll pitch
            out = obj.taitBryan312ToRotMat(rotVec312);
        end
        %
        function out = updateEuler321YawInRotMat(obj, yaw, rot_in)
            euler321 = obj.dcm2euler(rot_in);
            euler321(3) = yaw;
            out = obj.euler2dcm(euler321);
        end
        % 等价于matlab rotm2eul
        function out = dcm2euler(obj,dcm)
            theta = asin(-dcm(3, 1));
            
            if abs(theta - pi/2) < 1.0e-3
                phi = 0;
                psi = atan2(dcm(2, 3), dcm(1, 3));
                
            elseif abs(theta + pi/2) < 1.0e-3
                phi = 0;
                psi = atan2(-dcm(2, 3), -dcm(1, 3));
            else
                phi = atan2(dcm(3, 2), dcm(3, 3));
                psi = atan2(dcm(2, 1), dcm(1, 1));
            end
            out = [phi,theta,psi];
        end
        % 等价matlab eul2rotm，注意输入的euler顺序相反，matlab【psi-theta-phi】, euler2dcm
        %【phi-theta-psi】
        function dcm = euler2dcm(obj,euler)
            phi = euler(1);
            theta = euler(2);
            psi = euler(3);
            cosPhi = cos(phi);
            sinPhi = sin(phi);
            cosThe = cos(theta);
            sinThe = sin(theta);
            cosPsi = cos(psi);
            sinPsi = sin(psi);
            dcm = zeros(3,3);
            dcm(1, 1) = cosThe * cosPsi;
            dcm(1, 2) = -cosPhi * sinPsi + sinPhi * sinThe * cosPsi;
            dcm(1, 3) = sinPhi * sinPsi + cosPhi * sinThe * cosPsi;
            
            dcm(2, 1) = cosThe * sinPsi;
            dcm(2, 2) = cosPhi * cosPsi + sinPhi * sinThe * sinPsi;
            dcm(2, 3) = -sinPhi * cosPsi + cosPhi * sinThe * sinPsi;
            
            dcm(3, 1) = -sinThe;
            dcm(3, 2) = sinPhi * cosThe;
            dcm(3, 3) = cosPhi * cosThe;
        end
        % converts Tait-Bryan 312 sequence of rotations from frame 1 to frame 2
        % to the corresponding rotation matrix that rotates from frame 2 to frame 1
        % rot312(0) - First rotation is a RH rotation about the Z axis (rad)
        % rot312(1) - Second rotation is a RH rotation about the X axis (rad)
        % rot312(2) - Third rotation is a RH rotation about the Y axis (rad)
        % See http://www.atacolorado.com/eulersequences.doc
        function R = taitBryan312ToRotMat(obj, rot312)
            % Calculate the frame2 to frame 1 rotation matrix from a 312 Tait-Bryan rotation sequence
            c2 = cos(rot312(3)); % third rotation is pitch
            s2 = sin(rot312(3));
            s1 = sin(rot312(2)); % second rotation is roll
            c1 = cos(rot312(2));
            s0 = sin(rot312(1)); % first rotation is yaw
            c0 = cos(rot312(1));
            
            R = eye(3,3);
            R(1, 1) = c0 * c2 - s0 * s1 * s2;
            R(2, 2) = c0 * c1;
            R(3, 3) = c2 * c1;
            R(1, 2) = -c1 * s0;
            R(1, 3) = s2 * c0 + c2 * s1 * s0;
            R(2, 1) = c2 * s0 + s2 * s1 * c0;
            R(2, 3) = s0 * s2 - s1 * c0 * c2;
            R(3, 1) = -s2 * c1;
            R(3, 2) = s1;
        end
        
        % This function was autogenerated from a symbolic function. Do not modify by hand.
        % Symbolic function: yaw_est_predict_covariance
        % Args:
        %     state: Matrix31
        %     P: Matrix33
        %     d_vel: Matrix21
        %     d_vel_var: Scalar
        %     d_ang_var: Scalar
        % Outputs:
        %     P_new: Matrix33
        function P_new = YawEstPredictCovariance(obj, state, P, d_vel, d_vel_var,d_ang, d_ang_var)
            tmp0 = cos(state(3, 1));
            tmp1 = sin(state(3, 1));
            tmp2 = -tmp0 * d_vel(2, 1) - tmp1 * d_vel(1, 1);
            tmp3 = P(1, 3) + P(3, 3) * tmp2;
            tmp4 = tmp0^2 * d_vel_var + tmp1^2 * d_vel_var;
            tmp5 = tmp0 * d_vel(1, 1) - tmp1 * d_vel(2, 1);
            tmp6 = P(2, 3) + P(3, 3) * tmp5;
            tmp7 = d_ang^2 + 1;
            
            % Output terms (1)
            P_new = zeros(3);
            
            P_new(1, 1) = P(1, 1) + P(3, 1) * tmp2 + tmp2 * tmp3 + tmp4;
            P_new(2, 1) = 0;
            P_new(3, 1) = 0;
            P_new(1, 2) = P(1, 2) + P(3, 2) * tmp2 + tmp3 * tmp5;
            P_new(2, 2) = P(2, 2) + P(3, 2) * tmp5 + tmp4 + tmp5 * tmp6;
            P_new(3, 2) = 0;
            P_new(1, 3) = tmp3 * tmp7;
            P_new(2, 3) = tmp6 * tmp7;
            P_new(3, 3) = P(3, 3)* tmp7^2+d_ang_var;


        end
        % This function was autogenerated from a symbolic function. Do not modify by hand.
        % Symbolic function: yaw_est_compute_measurement_update
        % Args:
        %     P: Matrix33
        %     vel_obs_var: Scalar
        %     epsilon: Scalar
        % Outputs:
        %     S_inv: Matrix22
        %     S_det_inv: Scalar
        %     K: Matrix32
        %     P_new: Matrix33
        function [obj, S_inv, S_det_inv, K, P_new] = YawEstComputeMeasurementUpdate(obj, P, vel_obs_var, epsilon)
            % Total ops: 60
            
            % Input arrays
            
            % Intermediate terms (15)
            tmp0 = P(2, 2) + vel_obs_var;
            tmp1 = P(1, 1) + vel_obs_var;
            tmp2 = -P(1, 2) * P(2, 1) + tmp0 * tmp1;
            tmp3 = 1 / (tmp2 + epsilon * (2 * min(0, (((tmp2) > 0) - ((tmp2) < 0))) + 1));
            tmp4 = tmp0 * tmp3;
            tmp5 = P(2, 1) * tmp3;
            tmp6 = P(1, 2) * tmp3;
            tmp7 = tmp1 * tmp3;
            tmp8 = -P(1, 2) * tmp5;
            tmp9 = P(1, 1) * tmp4 + tmp8;
            tmp10 = -P(2, 2) * tmp5 + tmp0 * tmp5;
            tmp11 = P(3, 1) * tmp4 - P(3, 2) * tmp5;
            tmp12 = -P(1, 1) * tmp6 + tmp1 * tmp6;
            tmp13 = P(2, 2) * tmp7 + tmp8;
            tmp14 = -P(3, 1) * tmp6 + P(3, 2) * tmp7;
            
            % Output terms (4)
            S_inv = zeros(2,2);
            
            S_inv(1, 1) = tmp4;
            S_inv(2, 1) = -tmp5;
            S_inv(1, 2) = -tmp6;
            S_inv(2, 2) = tmp7;
            
            S_det_inv = tmp3;
            
            K = zeros(3,2);
            
            K(1, 1) = tmp9;
            K(2, 1) = tmp10;
            K(3, 1) = tmp11;
            K(1, 2) = tmp12;
            K(2, 2) = tmp13;
            K(3, 2) = tmp14;
            
            P_new = zeros(3,3);
            
            P_new(1, 1) = -P(1, 1) * tmp9 + P(1, 1) - P(2, 1) * tmp12;
            P_new(2, 1) = 0;
            P_new(3, 1) = 0;
            P_new(1, 2) = -P(1, 2) * tmp9 + P(1, 2) - P(2, 2) * tmp12;
            P_new(2, 2) = -P(1, 2) * tmp10 - P(2, 2) * tmp13 + P(2, 2);
            P_new(3, 2) = 0;
            P_new(1, 3) = -P(1, 3) * tmp9 + P(1, 3) - P(2, 3) * tmp12;
            P_new(2, 3) = -P(1, 3) * tmp10 - P(2, 3) * tmp13 + P(2, 3);
            P_new(3, 3) = -P(1, 3) * tmp11 - P(2, 3) * tmp14 + P(3, 3);
        end
        
 

        function  out = get_n_ekf_gsf(obj)
            out = obj.m_ekf_gsf;
        end

        function out = get_n_ahrs_ekf(obj)
            out = obj.m_ahrs_ekf_gsf;
        end


        function out = ger_n_weight(obj)
            out = obj.m_model_weights;
        end

        function gain = get_ahrs_fusion_gain(obj)
            gain = obj.m_ahrs_accel_fusion_gain;
        end
    end
end