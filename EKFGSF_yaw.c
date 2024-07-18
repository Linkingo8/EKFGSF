#include "EKFGSF_yaw.h"

// Parameters - these could be made tuneable
#define _gyro_noise     0.1f   // yaw rate noise used for covariance prediction (rad/sec)
#define _accel_noise    2.0f   // horizontal accel noise used for covariance prediction (m/sec**2)
#define _tilt_gain      0.2f   // gain from tilt error to gyro correction for complementary filter (1/sec)
#define _gyro_bias_gain 0.04f  // gain applied to integral of gyro correction for complementary filter (1/sec)
#define OFFSET          10

EKFGSF_YAW_t *_ekfgsf_yaw;

int nan_or_inf_identifier = NO_NAN_OR_INF_FOUND;

float nan_debug_data[16];

// calculate the gain from gravity vector misalingment to tilt correction to be used by all AHRS filters
static float ahrsCalcAccelGain(void);

// update specified AHRS rotation matrix using IMU and optionally true airspeed data
static void ahrsPredict(const uint8_t model_index, const Vector3f delta_ang, const float delta_ang_dt);

// align all AHRS roll and pitch orientations using IMU delta velocity vector
static void ahrsAlignTilt(const Vector3f delta_vel);

// align all AHRS yaw orientations to initial values
static void ahrsAlignYaw(void);

// Efficient propagation of a delta angle in body frame applied to the body to earth frame rotation matrix
static void ahrsPredictRotMat(e_matrix_t *R, const Vector3f g);

// initialise states and covariance data for the GSF and EKF filters
static void initialiseEKFGSF(const Vector2f vel_NE, const float vel_accuracy);

// predict state and covariance for the specified EKF using inertial data
static void predictEKF(const uint8_t model_index, const Vector3f delta_ang, const float delta_ang_dt, const Vector3f delta_vel, const float delta_vel_dt, bool in_air);

// update state and covariance for the specified EKF using a NE velocity measurement
// return false if update failed
static bool updateEKF(const uint8_t model_index, const Vector2f vel_NE, const float vel_accuracy);

// return the probability of the state estimate for the specified EKF assuming a gaussian error distribution
static float gaussianDensity(const uint8_t model_index);

static void YawEstPredictCovariance(const e_matrix_t *state,
                                    const e_matrix_t *P,
                                    const Vector2f    d_vel,
                                    const float       d_vel_var,
                                    const float       d_ang,
                                    const float       d_ang_var,
                                    e_matrix_t *const P_new);

static void YawEstComputeMeasurementUpdate(const e_matrix_t *P, const float vel_obs_var, const float epsilon, e_matrix_t *const S_inv, float *const S_det_inv, e_matrix_t *const K, e_matrix_t *const P_new);

static int checkIfNanOrInfInEKFGSF_YAW(const EKFGSF_YAW_t *ekfgsf_yaw, const unsigned int mark);

void EKFGSF_yaw(void)
{
    _ekfgsf_yaw = (EKFGSF_YAW_t *)pvPortMalloc(sizeof(EKFGSF_YAW_t));

    Vector3f zero_vector3f = {0.0f, 0.0f, 0.0f};
    // 需确保_ekfgsf_yaw所有成员都被正确初始化
    for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
    {
        _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R = (e_matrix_t *)pvPortMalloc(sizeof(e_matrix_t));
        InitMatrix(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R, 3, 3, 0.0f);
        _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias = zero_vector3f;

        _ekfgsf_yaw->_ekf_gsf[model_index].X         = (e_matrix_t *)pvPortMalloc(sizeof(e_matrix_t));
        _ekfgsf_yaw->_ekf_gsf[model_index].P         = (e_matrix_t *)pvPortMalloc(sizeof(e_matrix_t));
        _ekfgsf_yaw->_ekf_gsf[model_index].S_inverse = (e_matrix_t *)pvPortMalloc(sizeof(e_matrix_t));
        _ekfgsf_yaw->_ekf_gsf[model_index].innov     = (e_matrix_t *)pvPortMalloc(sizeof(e_matrix_t));

        InitMatrix(_ekfgsf_yaw->_ekf_gsf[model_index].X, 3, 1, 0.0f);
        InitMatrix(_ekfgsf_yaw->_ekf_gsf[model_index].P, 3, 3, 0.0f);
        InitMatrix(_ekfgsf_yaw->_ekf_gsf[model_index].S_inverse, 2, 2, 0.0f);
        InitMatrix(_ekfgsf_yaw->_ekf_gsf[model_index].innov, 2, 1, 0.0f);
        _ekfgsf_yaw->_ekf_gsf[model_index].S_det_inverse = 0.0f;

        _ekfgsf_yaw->_model_weights[model_index] = 0.0f;
    }

    _ekfgsf_yaw->_ahrs_accel                = zero_vector3f;
    _ekfgsf_yaw->_ahrs_ekf_gsf_tilt_aligned = false;
    _ekfgsf_yaw->_ekf_gsf_vel_fuse_started  = false;
    _ekfgsf_yaw->_gsf_yaw                   = 0.0f;
    _ekfgsf_yaw->_gsf_yaw_variance          = Sq(M_PI_F);
    _ekfgsf_yaw->_true_airspeed             = 0.0f;
}

void EKFGSF_yaw_predict(const imuSample_t imu_sample, bool in_air)
{
    if (nan_or_inf_identifier != 0)
        return;
    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    // 输入数据保护
    if (isnan(imu_sample.delta_ang.x) || isnan(imu_sample.delta_ang.y) || isnan(imu_sample.delta_ang.z) ||
        isnan(imu_sample.delta_vel.x) || isnan(imu_sample.delta_vel.y) || isnan(imu_sample.delta_vel.z) ||
        isinf(imu_sample.delta_ang.x) || isinf(imu_sample.delta_ang.y) || isinf(imu_sample.delta_ang.z) ||
        isinf(imu_sample.delta_vel.x) || isinf(imu_sample.delta_vel.y) || isinf(imu_sample.delta_vel.z))
    {
        return;
    }
    const Vector3f accel = Vector3f_divied(imu_sample.delta_vel, delta_dt);

    if (delta_dt > 0.001f)
    {
        // to reduce effect of vibration, filter using an LPF whose time constant is 1/10 of the AHRS tilt correction time constant
        const float filter_coef  = fminf(10.f * delta_dt * _tilt_gain, 1.0f);
        _ekfgsf_yaw->_ahrs_accel = Vector3f_Add(Vector3f_Scalar(accel, (1.0f - filter_coef)), Vector3f_Scalar(accel, filter_coef));
    }
    else
    {
        return;
    }

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    // Initialise states first time
    if (!_ekfgsf_yaw->_ahrs_ekf_gsf_tilt_aligned)
    {
        // check for excessive acceleration to reduce likelihood of large initial roll/pitch errors
        // due to vehicle movement
        const float accel_norm_sq     = Sq(Pythagorous3(accel.x, accel.y, accel.z));
        const float accel_lpf_norm_sq = Sq(Pythagorous3(_ekfgsf_yaw->_ahrs_accel.x, _ekfgsf_yaw->_ahrs_accel.y, _ekfgsf_yaw->_ahrs_accel.z));

        static const float upper_accel_limit = ONE_G * 1.1f;
        static const float lower_accel_limit = ONE_G * 0.9f;

        const bool ok_to_align = (accel_norm_sq > Sq(lower_accel_limit)) && (accel_norm_sq < Sq(upper_accel_limit)) && (accel_lpf_norm_sq > Sq(lower_accel_limit)) && (accel_lpf_norm_sq < Sq(upper_accel_limit));

        if (ok_to_align)
        {
            ahrsAlignTilt(imu_sample.delta_vel);
            _ekfgsf_yaw->_ahrs_ekf_gsf_tilt_aligned = true;

            nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
            if (nan_or_inf_identifier != 0)
                return;
        }
        else
        {
            return;
        }
    }

    for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
    {
        predictEKF(model_index, imu_sample.delta_ang, delta_dt, imu_sample.delta_vel, delta_dt, in_air);

        if (nan_or_inf_identifier != 0)
            return;
        nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
        if (nan_or_inf_identifier != 0)
            return;
    }
}

void EKFGSF_yaw_fuseVelocity(const Vector2f vel_NE, const float vel_accuracy, const bool in_air)
{
    if (nan_or_inf_identifier != 0)
        return;
    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    if (isnan(vel_NE.x) || isnan(vel_NE.y) || isnan(vel_accuracy) || isinf(vel_NE.x) || isinf(vel_NE.y) || isinf(vel_accuracy))
    {
        return;
    }

    // we don't start running the EKF part of the algorithm until there are regular velocity observations
    if (!_ekfgsf_yaw->_ekf_gsf_vel_fuse_started)
    {

        initialiseEKFGSF(vel_NE, vel_accuracy);

        nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
        if (nan_or_inf_identifier != 0)
            return;

        ahrsAlignYaw();

        nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
        if (nan_or_inf_identifier != 0)
            return;

        // don't start until in air or velocity is not negligible
        if (in_air || Pythagorous2(vel_NE.x, vel_NE.y) > vel_accuracy)
        {
            _ekfgsf_yaw->_ekf_gsf_vel_fuse_started = true;
        }
    }
    else
    {
        bool bad_update = false;

        for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
        {
            // subsequent measurements are fused as direct state observations
            if (!updateEKF(model_index, vel_NE, vel_accuracy))
            {
                bad_update = true;
            }
            nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
            if (nan_or_inf_identifier != 0)
                return;
        }

        if (!bad_update)
        {
            float total_weight = 0.0f;
            // calculate weighting for each model assuming a normal distribution
            const float min_weight     = 1e-5f;
            uint8_t     n_weight_clips = 0;

            for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
            {
                _ekfgsf_yaw->_model_weights[model_index] = gaussianDensity(model_index) * _ekfgsf_yaw->_model_weights[model_index];

                nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
                if (nan_or_inf_identifier != 0)
                    return;

                if (_ekfgsf_yaw->_model_weights[model_index] < min_weight)
                {
                    n_weight_clips++;
                    _ekfgsf_yaw->_model_weights[model_index] = min_weight;
                }

                nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
                if (nan_or_inf_identifier != 0)
                    return;

                total_weight += _ekfgsf_yaw->_model_weights[model_index];
            }

            // normalise the weighting function
            if (n_weight_clips < N_MODELS_EKFGSF)
            {
                for (uint8_t i = 0; i < N_MODELS_EKFGSF; i++)
                {
                    _ekfgsf_yaw->_model_weights[i] /= total_weight;

                    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
                    if (nan_or_inf_identifier != 0)
                        return;
                }
            }
            else
            {
                // all weights have collapsed due to excessive innovation variances so reset filters
                EKFGSF_yaw_reset();
            }
        }

        // Calculate a composite yaw vector as a weighted average of the states for each model.
        // To avoid issues with angle wrapping, the yaw state is converted to a vector with length
        // equal to the weighting value before it is summed.
        Vector2f yaw_vector = {0.0f, 0.0f};

        nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
        if (nan_or_inf_identifier != 0)
            return;

        for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
        {
            yaw_vector.x += _ekfgsf_yaw->_model_weights[model_index] * cosf(_ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0]);
            yaw_vector.y += _ekfgsf_yaw->_model_weights[model_index] * sinf(_ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0]);
        }

        _ekfgsf_yaw->_gsf_yaw = atan2f(yaw_vector.y, yaw_vector.x);

        nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
        if (nan_or_inf_identifier != 0)
            return;

        // calculate a composite variance for the yaw state from a weighted average of the variance for each model
        // models with larger innovations are weighted less
        _ekfgsf_yaw->_gsf_yaw_variance = 0.0f;

        for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
        {
            const float yaw_delta = WrapRadiansPi(_ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0] - _ekfgsf_yaw->_gsf_yaw);
            _ekfgsf_yaw->_gsf_yaw_variance += _ekfgsf_yaw->_model_weights[model_index] * (_ekfgsf_yaw->_ekf_gsf[model_index].P->mat[2][2] + yaw_delta * yaw_delta);

            nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
            if (nan_or_inf_identifier != 0)
                return;
        }
    }
}

void EKFGSF_yaw_reset(void)
{
    _ekfgsf_yaw->_ekf_gsf_vel_fuse_started = false;

    _ekfgsf_yaw->_gsf_yaw_variance = Sq(M_PI_F);
}

static void ahrsPredict(const uint8_t model_index, const Vector3f delta_ang, const float delta_ang_dt)
{
    // generate attitude solution using simple complementary filter for the selected model
    const Vector3f ang_rate = Vector3f_Sub(Vector3f_divied(delta_ang, fmaxf(delta_ang_dt, 0.001f)),
                                           _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias);

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    e_matrix_t R_to_body;
    MatrixTranspose(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R, &R_to_body);

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    Vector3f gravity_direction_bf;

    gravity_direction_bf.x = R_to_body.mat[0][2];
    gravity_direction_bf.y = R_to_body.mat[1][2];
    gravity_direction_bf.z = R_to_body.mat[2][2];

    DeInitMatrix(&R_to_body);
    const float ahrs_accel_norm = fmaxf(Pythagorous3(_ekfgsf_yaw->_ahrs_accel.x, _ekfgsf_yaw->_ahrs_accel.y, _ekfgsf_yaw->_ahrs_accel.z), 1e-5f);

    // gain from accel vector tilt error to rate gyro correction used by AHRS calculation
    const float ahrs_accel_fusion_gain = ahrsCalcAccelGain();
    nan_debug_data[6]                  = ahrs_accel_fusion_gain;
    // Perform angular rate correction using accel data and reduce correction as accel magnitude moves away from 1 g (reduces drift when vehicle picked up and moved).
    // During fixed wing flight, compensate for centripetal acceleration assuming coordinated turns and X axis forward
    Vector3f tilt_correction = {0.0f, 0.0f, 0.0f};

    if (ahrs_accel_fusion_gain > 0.f)
    {

        Vector3f accel = _ekfgsf_yaw->_ahrs_accel;

        if (isfinite(_ekfgsf_yaw->_true_airspeed) && (_ekfgsf_yaw->_true_airspeed > FLT_EPSILON))
        {
            // Calculate body frame centripetal acceleration with assumption X axis is aligned with the airspeed vector
            // Use cross product of body rate and body frame airspeed vector
            Vector3f centripetal_accel_bf;
            centripetal_accel_bf.x = 0.0f;
            centripetal_accel_bf.y = _ekfgsf_yaw->_true_airspeed * ang_rate.z;
            centripetal_accel_bf.z = -_ekfgsf_yaw->_true_airspeed * ang_rate.y;

            // correct measured accel for centripetal acceleration
            accel = Vector3f_Sub(accel, centripetal_accel_bf);
        }

        tilt_correction = Vector3f_Scalar(VectorCrossProduct(gravity_direction_bf, accel), ahrs_accel_fusion_gain / ahrs_accel_norm);
    }
    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    // Gyro bias estimation
    const float gyro_bias_limit = 0.05f;
    const float spin_rate       = Pythagorous3(ang_rate.x, ang_rate.y, ang_rate.z);

    if (spin_rate < 10.0f * DEG_TO_RAD)
    {
        _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias   = Vector3f_Sub(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias,
                                                                         Vector3f_Scalar(tilt_correction, (_gyro_bias_gain * delta_ang_dt)));
        _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias.x = ConstrainFloat(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias.x,
                                                                             -gyro_bias_limit,
                                                                             gyro_bias_limit);
        _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias.y = ConstrainFloat(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias.y,
                                                                             -gyro_bias_limit,
                                                                             gyro_bias_limit);
        _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias.z = ConstrainFloat(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias.z,
                                                                             -gyro_bias_limit,
                                                                             gyro_bias_limit);
    }
    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    // delta angle from previous to current frame
    const Vector3f delta_angle_corrected = Vector3f_Add(delta_ang,
                                                        Vector3f_Scalar(Vector3f_Sub(tilt_correction, _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias), delta_ang_dt));

    // Apply delta angle to rotation matrix
    ahrsPredictRotMat(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R, delta_angle_corrected);
    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;
}

static void ahrsAlignTilt(const Vector3f delta_vel)
{
    // Rotation matrix is constructed directly from acceleration measurement and will be the same for
    // all models so only need to calculate it once. Assumptions are:
    // 1) Yaw angle is zero - yaw is aligned later for each model when velocity fusion commences.
    // 2) The vehicle is not accelerating so all of the measured acceleration is due to gravity.

    // Calculate earth frame Down axis unit vector rotated into body frame
    Vector3f d_vel_ = delta_vel;
    Vector3f_Normalize(&d_vel_);
    const Vector3f down_in_bf = Vector3f_Scalar(d_vel_, -1);

    // Calculate earth frame North axis unit vector rotated into body frame, orthogonal to 'down_in_bf'
    Vector3f i_vec_bf;
    i_vec_bf.x = 1.0f;
    i_vec_bf.y = 0.0f;
    i_vec_bf.z = 0.0f;

    Vector3f north_in_bf = Vector3f_Sub(i_vec_bf, Vector3f_Scalar(down_in_bf, VectorDotProduct(i_vec_bf, down_in_bf)));

    Vector3f_Normalize(&north_in_bf);

    // Calculate earth frame East axis unit vector rotated into body frame, orthogonal to 'down_in_bf' and 'north_in_bf'
    const Vector3f east_in_bf = VectorCrossProduct(down_in_bf, north_in_bf);

    // Each column in a rotation matrix from earth frame to body frame represents the projection of the
    // corresponding earth frame unit vector rotated into the body frame, eg 'north_in_bf' would be the first column.
    // We need the rotation matrix from body frame to earth frame so the earth frame unit vectors rotated into body
    // frame are copied into corresponding rows instead.
    e_matrix_t R;
    InitMatrix(&R, 3, 3, 0.0f);
    R.mat[0][0] = north_in_bf.x;
    R.mat[0][1] = north_in_bf.y;
    R.mat[0][2] = north_in_bf.z;
    R.mat[1][0] = east_in_bf.x;
    R.mat[1][1] = east_in_bf.y;
    R.mat[1][2] = east_in_bf.z;
    R.mat[2][0] = down_in_bf.x;
    R.mat[2][1] = down_in_bf.y;
    R.mat[2][2] = down_in_bf.z;

    for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
    {
        CopyMatrixValues(&R, _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R);
    }
    DeInitMatrix(&R);
}

static void ahrsAlignYaw(void)
{
    // Align yaw angle for each model
    for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
    {

        const float yaw = WrapRadiansPi(_ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0]);
        updateYawInRotMat(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R, yaw);

        nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
        if (nan_or_inf_identifier != 0)
            return;
    }
}

static void predictEKF(const uint8_t model_index, const Vector3f delta_ang, const float delta_ang_dt, const Vector3f delta_vel, const float delta_vel_dt, bool in_air)
{
    // generate an attitude reference using IMU data
    ahrsPredict(model_index, delta_ang, delta_ang_dt);

    if (nan_or_inf_identifier != 0)
        return;
    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    // we don't start running the EKF part of the algorithm until there are regular velocity observations
    if (!_ekfgsf_yaw->_ekf_gsf_vel_fuse_started)
    {
        return;
    }

    // Calculate the yaw state using a projection onto the horizontal that avoids gimbal lock
    _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0] = getEulerYaw(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R);

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    // calculate delta velocity in a horizontal front-right frame
    const Vector3f del_vel_NED = eMatrix3MulVector3(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R, delta_vel);
    const float    cos_yaw     = cosf(_ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0]);
    const float    sin_yaw     = sinf(_ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0]);
    const float    dvx         = del_vel_NED.x * cos_yaw + del_vel_NED.y * sin_yaw;
    const float    dvy         = -del_vel_NED.x * sin_yaw + del_vel_NED.y * cos_yaw;

    if (nan_or_inf_identifier != 0)
        return;

    Vector3f    del_ang_NED = eMatrix3MulVector3(_ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R, delta_ang);
    const float daz         = del_ang_NED.z;

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    // delta velocity process noise double if we're not in air
    const float accel_noise = in_air ? _accel_noise : 2.0f * _accel_noise;
    const float d_vel_var   = Sq(accel_noise * delta_vel_dt);

    // Use fixed values for delta angle process noise variances
    const float d_ang_var = Sq(_gyro_noise * delta_ang_dt);

    Vector2f vel_b;
    vel_b.x = dvx;
    vel_b.y = dvy;

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    YawEstPredictCovariance(_ekfgsf_yaw->_ekf_gsf[model_index].X, _ekfgsf_yaw->_ekf_gsf[model_index].P, vel_b, d_vel_var, daz, d_ang_var, _ekfgsf_yaw->_ekf_gsf[model_index].P);
    if (nan_or_inf_identifier != 0)
        return;

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;
    // covariance matrix is symmetrical, so copy upper half to lower half
    _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[1][0] = _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[0][1];
    _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[2][0] = _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[0][2];
    _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[2][1] = _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[1][2];

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;
    // constrain variances
    const float min_var = 1e-6f;

    for (unsigned index = 0; index < 3; index++)
    {
        _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[index][index] = fmaxf(_ekfgsf_yaw->_ekf_gsf[model_index].P->mat[index][index], min_var);
    }

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;

    // sum delta velocities in earth frame:
    _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[0][0] += del_vel_NED.x;
    _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[1][0] += del_vel_NED.y;

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return;
}

static bool updateEKF(const uint8_t model_index, const Vector2f vel_NE, const float vel_accuracy)
{
    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    // set observation variance from accuracy estimate supplied by GPS and apply a sanity check minimum
    const float vel_obs_var = Sq(fmaxf(vel_accuracy, 0.01f));

    // calculate velocity observation innovations
    _ekfgsf_yaw->_ekf_gsf[model_index].innov->mat[0][0] = _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[0][0] - vel_NE.x;
    _ekfgsf_yaw->_ekf_gsf[model_index].innov->mat[1][0] = _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[1][0] - vel_NE.y;

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    e_matrix_t K;

    InitMatrix(&K, 3, 2, 0.0f);

    YawEstComputeMeasurementUpdate(_ekfgsf_yaw->_ekf_gsf[model_index].P,
                                   vel_obs_var,
                                   FLT_EPSILON,
                                   _ekfgsf_yaw->_ekf_gsf[model_index].S_inverse,
                                   &_ekfgsf_yaw->_ekf_gsf[model_index].S_det_inverse,
                                   &K,
                                   _ekfgsf_yaw->_ekf_gsf[model_index].P);

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    // copy upper to lower diagonal
    _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[1][0] = _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[0][1];
    _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[2][0] = _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[0][2];
    _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[2][1] = _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[1][2];

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    // constrain variances
    const float min_var = 1e-6f;

    for (unsigned index = 0; index < 3; index++)
    {
        _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[index][index] = fmaxf(_ekfgsf_yaw->_ekf_gsf[model_index].P->mat[index][index], min_var);
    }

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    // test ratio = transpose(innovation) * inverse(innovation variance) * innovation = [1x2] * [2,2] * [2,1] = [1,1]
    e_matrix_t innov_transpose;
    MatrixTranspose(_ekfgsf_yaw->_ekf_gsf[model_index].innov, &innov_transpose);
    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    e_matrix_t tmp1, res;
    MatrixMultiplication(_ekfgsf_yaw->_ekf_gsf[model_index].S_inverse, _ekfgsf_yaw->_ekf_gsf[model_index].innov, &tmp1);

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    MatrixMultiplication(&innov_transpose, &tmp1, &res);

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    const float test_ratio = res.mat[0][0];

    DeInitMatrix(&innov_transpose);
    DeInitMatrix(&tmp1);
    DeInitMatrix(&res);

    // Perform a chi-square innovation consistency test and calculate a compression scale factor
    // that limits the magnitude of innovations to 5-sigma
    // If the test ratio is greater than 25 (5 Sigma) then reduce the length of the innovation vector to clip it at 5-Sigma
    // This protects from large measurement spikes
    const float innov_comp_scale_factor = test_ratio > 25.0f ? sqrtf(25.0f / test_ratio) : 1.0f;

    // Correct the state vector and capture the change in yaw angle
    const float oldYaw = _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0];

    e_matrix_t tmp2, tmp3, tmp4;
    MatrixMultiplication(&K, _ekfgsf_yaw->_ekf_gsf[model_index].innov, &tmp2);

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    MatrixScalarMultiplication(&tmp2, innov_comp_scale_factor, &tmp3);

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    MatrixSubtraction(_ekfgsf_yaw->_ekf_gsf[model_index].X, &tmp3, &tmp4);

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    CopyMatrixValues(&tmp4, _ekfgsf_yaw->_ekf_gsf[model_index].X);

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    DeInitMatrix(&tmp2);
    DeInitMatrix(&tmp3);
    DeInitMatrix(&tmp4);

    const float yawDelta = _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0] - oldYaw;

    // apply the change in yaw angle to the AHRS
    // take advantage of sparseness in the yaw rotation matrix
    const float cosYaw   = cosf(yawDelta);
    const float sinYaw   = sinf(yawDelta);
    const float R_prev00 = _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[0][0];
    const float R_prev01 = _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[0][1];
    const float R_prev02 = _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[0][2];

    _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[0][0] = R_prev00 * cosYaw - _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[1][0] * sinYaw;
    _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[0][1] = R_prev01 * cosYaw - _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[1][1] * sinYaw;
    _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[0][2] = R_prev02 * cosYaw - _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[1][2] * sinYaw;
    _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[1][0] = R_prev00 * sinYaw + _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[1][0] * cosYaw;
    _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[1][1] = R_prev01 * sinYaw + _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[1][1] * cosYaw;
    _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[1][2] = R_prev02 * sinYaw + _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].R->mat[1][2] * cosYaw;

    nan_or_inf_identifier = checkIfNanOrInfInEKFGSF_YAW(_ekfgsf_yaw, __LINE__);
    if (nan_or_inf_identifier != 0)
        return false;

    DeInitMatrix(&K);
    return true;
}

static void initialiseEKFGSF(const Vector2f vel_NE, const float vel_accuracy)
{
    _ekfgsf_yaw->_gsf_yaw          = 0.0f;
    _ekfgsf_yaw->_gsf_yaw_variance = Sq(M_PI_F / 2.0f);

    for (uint8_t i = 0; i < N_MODELS_EKFGSF; i++)
    {
        _ekfgsf_yaw->_model_weights[i] = (1.0f / (float)N_MODELS_EKFGSF);  // All filter models start with the same weight
    }

    const float yaw_increment = 2.0f * M_PI_F / (float)N_MODELS_EKFGSF;

    for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
    {

        _ekfgsf_yaw->_ekf_gsf[model_index].S_det_inverse = 0.0f;

        // evenly space initial yaw estimates in the region between +-Pi
        _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0] = -M_PI_F + (0.5f * yaw_increment) + ((float)model_index * yaw_increment);

        // take velocity states and corresponding variance from last measurement
        _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[0][0] = vel_NE.x;
        _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[1][0] = vel_NE.y;

        _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[0][0] = Sq(fmaxf(vel_accuracy, 0.01f));
        _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[1][1] = Sq(fmaxf(vel_accuracy, 0.01f));

        // use half yaw interval for yaw uncertainty
        _ekfgsf_yaw->_ekf_gsf[model_index].P->mat[2][2] = Sq(0.5f * yaw_increment);
    }
}

static float gaussianDensity(const uint8_t model_index)
{
    // calculate transpose(innovation) * inv(S) * innovation
    e_matrix_t innov_transpose;
    MatrixTranspose(_ekfgsf_yaw->_ekf_gsf[model_index].innov, &innov_transpose);

    e_matrix_t tmp_1, res;
    MatrixMultiplication(_ekfgsf_yaw->_ekf_gsf[model_index].S_inverse, _ekfgsf_yaw->_ekf_gsf[model_index].innov, &tmp_1);
    MatrixMultiplication(&innov_transpose, &tmp_1, &res);

    const float normDist = res.mat[0][0];

    DeInitMatrix(&innov_transpose);
    DeInitMatrix(&tmp_1);
    DeInitMatrix(&res);

    return (1.0f / (2.0f * M_PI_F)) * sqrtf(_ekfgsf_yaw->_ekf_gsf[model_index].S_det_inverse) * expf(-0.5f * normDist);
}

bool EKFGSF_yaw_getLogData(float *yaw_composite, float *yaw_variance, float yaw[N_MODELS_EKFGSF], float innov_VN[N_MODELS_EKFGSF], float innov_VE[N_MODELS_EKFGSF], float weight[N_MODELS_EKFGSF], float debug_data[16])
{
    // if (_ekfgsf_yaw->_ekf_gsf_vel_fuse_started)
    {
        *yaw_composite = _ekfgsf_yaw->_gsf_yaw;
        *yaw_variance  = _ekfgsf_yaw->_gsf_yaw_variance;

        for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
        {
            yaw[model_index]      = _ekfgsf_yaw->_ekf_gsf[model_index].X->mat[2][0];
            innov_VN[model_index] = _ekfgsf_yaw->_ekf_gsf[model_index].innov->mat[0][0];
            innov_VE[model_index] = _ekfgsf_yaw->_ekf_gsf[model_index].innov->mat[1][0];
            weight[model_index]   = _ekfgsf_yaw->_model_weights[model_index];
        }

        for (uint8_t k = 0; k < 16; k++)
        {
            debug_data[k] = nan_debug_data[k];
        }

        return true;
    }
}

static float ahrsCalcAccelGain(void)
{
    // Calculate the acceleration fusion gain using a continuous function that is unity at 1g and zero
    // at the min and max g value. Allow for more acceleration when flying as a fixed wing vehicle using centripetal
    // acceleration correction as higher and more sustained g will be experienced.
    // Use a quadratic instead of linear function to prevent vibration around 1g reducing the tilt correction effectiveness.
    // see https://www.desmos.com/calculator/dbqbxvnwfg

    float      attenuation                            = 2.0f;
    const bool centripetal_accel_compensation_enabled = isfinite(_ekfgsf_yaw->_true_airspeed) && (_ekfgsf_yaw->_true_airspeed > FLT_EPSILON);

    const float ahrs_accel_norm = Pythagorous3(_ekfgsf_yaw->_ahrs_accel.x, _ekfgsf_yaw->_ahrs_accel.y, _ekfgsf_yaw->_ahrs_accel.z);

    if (centripetal_accel_compensation_enabled && (ahrs_accel_norm > ONE_G))
    {
        attenuation = 1.0f;
    }

    const float delta_accel_g = (ahrs_accel_norm - ONE_G) / ONE_G;
    return _tilt_gain * Sq(1.0f - min(attenuation * fabsf(delta_accel_g), 1.0f));
}

static void ahrsPredictRotMat(e_matrix_t *R, const Vector3f g)
{
    e_matrix_t R_o;
    InitMatrix(&R_o, 3, 3, 0.0f);
    CopyMatrixValues(R, &R_o);

    R->mat[0][0] += R_o.mat[0][1] * g.z - R_o.mat[0][2] * g.y;
    R->mat[0][1] += R_o.mat[0][2] * g.x - R_o.mat[0][0] * g.z;
    R->mat[0][2] += R_o.mat[0][0] * g.y - R_o.mat[0][1] * g.x;
    R->mat[1][0] += R_o.mat[1][1] * g.z - R_o.mat[1][2] * g.y;
    R->mat[1][1] += R_o.mat[1][2] * g.x - R_o.mat[1][0] * g.z;
    R->mat[1][2] += R_o.mat[1][0] * g.y - R_o.mat[1][1] * g.x;
    R->mat[2][0] += R_o.mat[2][1] * g.z - R_o.mat[2][2] * g.y;
    R->mat[2][1] += R_o.mat[2][2] * g.x - R_o.mat[2][0] * g.z;
    R->mat[2][2] += R_o.mat[2][0] * g.y - R_o.mat[2][1] * g.x;

    // Renormalise rows
    for (uint8_t r = 0; r < 3; r++)
    {
        const float rowLengthSq = Pythagorous3(R->mat[r][0], R->mat[r][1], R->mat[r][2]);

        if (rowLengthSq > FLT_EPSILON)
        {
            // Use linear approximation for inverse sqrt taking advantage of the row length being close to 1.0
            const float rowLengthInv = 1.5f - 0.5f * rowLengthSq;
            R->mat[r][0] *= rowLengthInv;
            R->mat[r][1] *= rowLengthInv;
            R->mat[r][2] *= rowLengthInv;
        }
    }
    DeInitMatrix(&R_o);
}

bool EKFGSF_yaw_isActive(void)
{
    bool active = (_ekfgsf_yaw->_ekf_gsf_vel_fuse_started && (nan_or_inf_identifier == NO_NAN_OR_INF_FOUND));
    return active;
}

float EKFGSF_yaw_getYaw(void)
{
    float yaw;
    if (_ekfgsf_yaw->_ekf_gsf_vel_fuse_started)
    {
        yaw = _ekfgsf_yaw->_gsf_yaw;
    }
    else
    {
        yaw = 0.0f;
    }
    return yaw;
}

float EKFGSF_yaw_getYawVar(void)
{
    float variance;
    if (_ekfgsf_yaw->_ekf_gsf_vel_fuse_started)
    {
        variance = _ekfgsf_yaw->_gsf_yaw_variance;
    }
    else
    {
        variance = Sq(2 * M_PI);
    }
    return variance;
}

int EKFGSF_yaw_getNanCheckIdentifier(void)
{
    return nan_or_inf_identifier;
}

bool EKFGSF_yaw_getAhrsTiltAligned(void)
{
    return _ekfgsf_yaw->_ahrs_ekf_gsf_tilt_aligned;
}

void EKFGSF_yaw_setTrueAirspeed(float true_airspeed)
{
    _ekfgsf_yaw->_true_airspeed = true_airspeed;
}

void EKFGSF_yaw_setGyroBias(const Vector3f imu_gyro_bias)
{
    // Initialise to gyro bias estimate from main filter because there could be a large
    // uncorrected rate gyro bias error about the gravity vector
    if (!_ekfgsf_yaw->_ahrs_ekf_gsf_tilt_aligned || !_ekfgsf_yaw->_ekf_gsf_vel_fuse_started)
    {
        // init gyro bias for each model
        for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++)
        {
            _ekfgsf_yaw->_ahrs_ekf_gsf[model_index].gyro_bias = imu_gyro_bias;
        }
    }
}

static void YawEstPredictCovariance(const e_matrix_t *state,
                                    const e_matrix_t *P,
                                    const Vector2f    d_vel,
                                    const float       d_vel_var,
                                    const float       d_ang,
                                    const float       d_ang_var,
                                    e_matrix_t *const P_new)
{

    e_matrix_t p_old_tmp;
    InitMatrix(&p_old_tmp, 3, 3, 0.0f);
    CopyMatrixValues(P, &p_old_tmp);

    const float _tmp0 = cosf(state->mat[2][0]);
    const float _tmp1 = sinf(state->mat[2][0]);
    const float _tmp2 = -_tmp0 * d_vel.y - _tmp1 * d_vel.x;
    const float _tmp3 = p_old_tmp.mat[0][2] + p_old_tmp.mat[2][2] * _tmp2;
    const float _tmp4 = Sq(_tmp0) * d_vel_var + Sq(_tmp1) * d_vel_var;
    const float _tmp5 = _tmp0 * d_vel.x - _tmp1 * d_vel.y;
    const float _tmp6 = p_old_tmp.mat[1][2] + p_old_tmp.mat[2][2] * _tmp5;
    const float _tmp7 = Sq(d_ang) + 1;

    if (P_new != NULL)
    {
        P_new->mat[0][0] = p_old_tmp.mat[0][0] + p_old_tmp.mat[2][0] * _tmp2 + _tmp2 * _tmp3 + _tmp4;
        P_new->mat[1][0] = 0.0f;
        P_new->mat[2][0] = 0.0f;
        P_new->mat[0][1] = p_old_tmp.mat[0][1] + p_old_tmp.mat[2][1] * _tmp2 + _tmp3 * _tmp5;
        P_new->mat[1][1] = p_old_tmp.mat[1][1] + p_old_tmp.mat[2][1] * _tmp5 + _tmp4 + _tmp5 * _tmp6;
        P_new->mat[2][1] = 0.0f;
        P_new->mat[0][2] = _tmp3 * _tmp7;
        P_new->mat[1][2] = _tmp6 * _tmp7;
        P_new->mat[2][2] = p_old_tmp.mat[2][2] * Sq(_tmp7) + d_ang_var;
    }

    DeInitMatrix(&p_old_tmp);
}

static void YawEstComputeMeasurementUpdate(const e_matrix_t *P, const float vel_obs_var, const float epsilon, e_matrix_t *const S_inv, float *const S_det_inv, e_matrix_t *const K, e_matrix_t *const P_new)
{

    e_matrix_t p_old_tmp;
    InitMatrix(&p_old_tmp, 3, 3, 0.0f);
    CopyMatrixValues(P, &p_old_tmp);

    const float _tmp0  = p_old_tmp.mat[1][1] + vel_obs_var;
    const float _tmp1  = p_old_tmp.mat[0][0] + vel_obs_var;
    const float _tmp2  = -p_old_tmp.mat[0][1] * p_old_tmp.mat[1][0] + _tmp0 * _tmp1;
    const float _tmp3  = 1.0f / (_tmp2 + epsilon * (2 * min(0, (((_tmp2) > 0) - ((_tmp2) < 0))) + 1));
    const float _tmp4  = _tmp0 * _tmp3;
    const float _tmp5  = p_old_tmp.mat[1][0] * _tmp3;
    const float _tmp6  = p_old_tmp.mat[0][1] * _tmp3;
    const float _tmp7  = _tmp1 * _tmp3;
    const float _tmp8  = -p_old_tmp.mat[0][1] * _tmp5;
    const float _tmp9  = p_old_tmp.mat[0][0] * _tmp4 + _tmp8;
    const float _tmp10 = -p_old_tmp.mat[1][1] * _tmp5 + _tmp0 * _tmp5;
    const float _tmp11 = p_old_tmp.mat[2][0] * _tmp4 - p_old_tmp.mat[2][1] * _tmp5;
    const float _tmp12 = -p_old_tmp.mat[0][0] * _tmp6 + _tmp1 * _tmp6;
    const float _tmp13 = p_old_tmp.mat[1][1] * _tmp7 + _tmp8;
    const float _tmp14 = -p_old_tmp.mat[2][0] * _tmp6 + p_old_tmp.mat[2][1] * _tmp7;

    if (S_inv != NULL)
    {
        S_inv->mat[0][0] = _tmp4;
        S_inv->mat[1][0] = -_tmp5;
        S_inv->mat[0][1] = -_tmp6;
        S_inv->mat[1][1] = _tmp7;
    }

    if (S_det_inv != NULL)
        *S_det_inv = _tmp3;

    if (K != NULL)
    {
        K->mat[0][0] = _tmp9;
        K->mat[1][0] = _tmp10;
        K->mat[2][0] = _tmp11;
        K->mat[0][1] = _tmp12;
        K->mat[1][1] = _tmp13;
        K->mat[2][1] = _tmp14;
    }

    if (P_new != NULL)
    {
        P_new->mat[0][0] = -p_old_tmp.mat[0][0] * _tmp9 + p_old_tmp.mat[0][0] - p_old_tmp.mat[1][0] * _tmp12;
        P_new->mat[1][0] = 0.0f;
        P_new->mat[2][0] = 0.0f;
        P_new->mat[0][1] = -p_old_tmp.mat[0][1] * _tmp9 + p_old_tmp.mat[0][1] - p_old_tmp.mat[1][1] * _tmp12;
        P_new->mat[1][1] = -p_old_tmp.mat[0][1] * _tmp10 - p_old_tmp.mat[1][1] * _tmp13 + p_old_tmp.mat[1][1];
        P_new->mat[2][1] = 0.0f;
        P_new->mat[0][2] = -p_old_tmp.mat[0][2] * _tmp9 + p_old_tmp.mat[0][2] - p_old_tmp.mat[1][2] * _tmp12;
        P_new->mat[1][2] = -p_old_tmp.mat[0][2] * _tmp10 - p_old_tmp.mat[1][2] * _tmp13 + p_old_tmp.mat[1][2];
        P_new->mat[2][2] = -p_old_tmp.mat[0][2] * _tmp11 - p_old_tmp.mat[1][2] * _tmp14 + p_old_tmp.mat[2][2];
    }

    DeInitMatrix(&p_old_tmp);
}

// 检查单个浮点数是否为NaN或Inf
static int checkFloatForNaNOrInf(float value, int nanIdentifier, int infIdentifier)
{
    if (isnan(value))
        return nanIdentifier;
    if (isinf(value))
        return infIdentifier;
    return NO_NAN_OR_INF_FOUND;
}

// 检查Vector3f是否包含NaN或Inf
static int checkVector3fForNaNOrInf(const Vector3f *vec, int baseNanIdentifier, int baseInfIdentifier)
{
    if (isnan(vec->x))
        return baseNanIdentifier + 1;
    if (isnan(vec->y))
        return baseNanIdentifier + 2;
    if (isnan(vec->z))
        return baseNanIdentifier + 3;
    if (isinf(vec->x))
        return baseInfIdentifier + 1;
    if (isinf(vec->y))
        return baseInfIdentifier + 2;
    if (isinf(vec->z))
        return baseInfIdentifier + 3;
    return NO_NAN_OR_INF_FOUND;
}

// 检查矩阵是否包含NaN或Inf
static int checkMatrixForNaNOrInf(const e_matrix_t *matrix, int baseNanIdentifier, int baseInfIdentifier)
{
    for (int i = 0; i < matrix->rows; ++i)
    {
        for (int j = 0; j < matrix->cols; ++j)
        {
            if (isnan(matrix->mat[i][j]))
            {
                return baseNanIdentifier + (i * matrix->cols + j) + 1;
            }
            if (isinf(matrix->mat[i][j]))
            {
                return baseInfIdentifier + (i * matrix->cols + j) + 1;
            }
        }
    }
    return NO_NAN_OR_INF_FOUND;
}

// 检查EKFGSF_YAW_t结构体是否包含NaN或Inf
static int checkIfNanOrInfInEKFGSF_YAW(const EKFGSF_YAW_t *ekfgsf_yaw, const unsigned int mark)
{
    int result;

    // 检查单个浮点数成员
    result = checkFloatForNaNOrInf(ekfgsf_yaw->_true_airspeed, NAN_TRUE_AIRSPEED, INF_TRUE_AIRSPEED);
    if (result != NO_NAN_OR_INF_FOUND)
        return result + mark * 10000;

    result = checkFloatForNaNOrInf(ekfgsf_yaw->_gsf_yaw, NAN_GSF_YAW, INF_GSF_YAW);
    if (result != NO_NAN_OR_INF_FOUND)
        return result + mark * 10000;

    result = checkFloatForNaNOrInf(ekfgsf_yaw->_gsf_yaw_variance, NAN_GSF_YAW_VARIANCE, INF_GSF_YAW_VARIANCE);
    if (result != NO_NAN_OR_INF_FOUND)
        return result + mark * 10000;

    // 循环检查Vector3f成员和矩阵
    for (int i = 0; i < N_MODELS_EKFGSF; ++i)
    {
        result = checkMatrixForNaNOrInf(ekfgsf_yaw->_ahrs_ekf_gsf[i].R, NAN_AHRS_EKF_GSF_R + i * OFFSET, INF_AHRS_EKF_GSF_R + i * OFFSET);
        if (result != NO_NAN_OR_INF_FOUND)
            return result + mark * 10000;

        result = checkVector3fForNaNOrInf(&ekfgsf_yaw->_ahrs_ekf_gsf[i].gyro_bias, NAN_AHRS_EKF_GSF_GYRO_BIAS_X + i * OFFSET, INF_AHRS_EKF_GSF_GYRO_BIAS_X + i * OFFSET);
        if (result != NO_NAN_OR_INF_FOUND)
            return result + mark * 10000;

        // 检查结构体内的e_matrix_t*成员
        result = checkMatrixForNaNOrInf(ekfgsf_yaw->_ekf_gsf[i].X, NAN_EKF_GSF_X + i * OFFSET, INF_EKF_GSF_X + i * OFFSET);
        if (result != NO_NAN_OR_INF_FOUND)
            return result + mark * 10000;

        result = checkMatrixForNaNOrInf(ekfgsf_yaw->_ekf_gsf[i].P, NAN_EKF_GSF_P + i * OFFSET, INF_EKF_GSF_P + i * OFFSET);
        if (result != NO_NAN_OR_INF_FOUND)
            return result + mark * 10000;

        result = checkMatrixForNaNOrInf(ekfgsf_yaw->_ekf_gsf[i].S_inverse, NAN_EKF_GSF_S_INVERSE + i * OFFSET, INF_EKF_GSF_S_INVERSE + i * OFFSET);
        if (result != NO_NAN_OR_INF_FOUND)
            return result + mark * 10000;

        // 检查结构体内的单个浮点数
        result = checkFloatForNaNOrInf(ekfgsf_yaw->_ekf_gsf[i].S_det_inverse, NAN_EKF_GSF_S_DET_INVERSE + i, INF_EKF_GSF_S_DET_INVERSE + i);
        if (result != NO_NAN_OR_INF_FOUND)
            return result + mark * 10000;
    }

    // 检查Vector3f _ahrs_accel
    result = checkVector3fForNaNOrInf(&ekfgsf_yaw->_ahrs_accel, NAN_AHRS_ACCEL_X, INF_AHRS_ACCEL_X);
    if (result != NO_NAN_OR_INF_FOUND)
        return result + mark * 10000;

    // 检查模型权重
    for (int i = 0; i < N_MODELS_EKFGSF; ++i)
    {
        result = checkFloatForNaNOrInf(ekfgsf_yaw->_model_weights[i], NAN_MODEL_WEIGHTS + i, INF_MODEL_WEIGHTS + i);
        if (result != NO_NAN_OR_INF_FOUND)
            return result + mark * 10000;
    }

    return NO_NAN_OR_INF_FOUND;  // 表明未找到NaN或Inf值
}
