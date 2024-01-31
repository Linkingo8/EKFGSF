#ifndef EKF_EKFGSF_YAW_H
#define EKF_EKFGSF_YAW_H

#include "vector3.h"
#include "mathTool.h"
#include "math.h"
#include "float.h"

#define  N_MODELS_EKFGSF  5
#define  ONE_G            9.80665f	
#define  delta_dt         0.004f

#define GSF_YAW_SAFE_ACCURACY_DEG 25.0f
#define GSF_YAW_GOOD_ACCURACY_DEG 15.0f

typedef struct {
	Vector3f    delta_ang;              ///< delta angle in body frame (integrated gyro measurements) (rad)
	Vector3f    delta_vel;              ///< delta velocity in body frame (integrated accelerometer measurements) (m/sec)
}imuSample_t;

//public:
	void EKFGSF_yaw(void);

	// Update Filter States - this should be called whenever new IMU data is available
	void EKFGSF_yaw_predict(const imuSample_t imu_sample, bool in_air);

	void EKFGSF_yaw_fuseVelocity(const Vector2f vel_NE, // NE velocity measurement (m/s)
			  const float vel_accuracy,	  // 1-sigma accuracy of velocity measurement (m/s)
			  bool in_air);

	void EKFGSF_yaw_setTrueAirspeed(float true_airspeed);
	void EKFGSF_yaw_setGyroBias(const Vector3f imu_gyro_bias);

	// get solution data for logging
	bool EKFGSF_yaw_getLogData(float *yaw_composite,
			                    float *yaw_composite_variance,
								float yaw[N_MODELS_EKFGSF],
								float innov_VN[N_MODELS_EKFGSF],
								float innov_VE[N_MODELS_EKFGSF],
								float weight[N_MODELS_EKFGSF]);

	bool  EKFGSF_yaw_isActive(void);
	float EKFGSF_yaw_getYaw(void);
	float EKFGSF_yaw_getYawVar(void);
	void  EKFGSF_yaw_reset(void);

///private
	typedef struct
	{
		// Declarations used by the bank of N_MODELS_EKFGSF AHRS complementary filters
		float _true_airspeed;				// true airspeed used for centripetal accel compensation (m/s)

		struct {
			e_matrix_t* R; 						// matrix that rotates a vector from body to earth frame
			Vector3f 	gyro_bias;            	// gyro bias learned and used by the quaternion calculation
		} _ahrs_ekf_gsf[N_MODELS_EKFGSF];

		bool _ahrs_ekf_gsf_tilt_aligned;  // true the initial tilt alignment has been calculated
		Vector3f _ahrs_accel;     // low pass filtered body frame specific force vector used by AHRS calculation (m/s/s)

		// Declarations used by a bank of N_MODELS_EKFGSF EKFs
		struct _ekf_gsf_struct {
			e_matrix_t* X;                       // Vel North (m/s),  Vel East (m/s), yaw (rad)s
			e_matrix_t* P;         				// covariance matrix
			e_matrix_t* S_inverse; 				// inverse of the innovation covariance matrix
			float S_det_inverse;                // inverse of the innovation covariance matrix determinant
			e_matrix_t* innov;                   // Velocity N,E innovation (m/s)
		} _ekf_gsf[N_MODELS_EKFGSF];

		bool _ekf_gsf_vel_fuse_started; // true when the EKF's have started fusing velocity data and the prediction and update processing is active

		// Declarations used by the Gaussian Sum Filter (GSF) that combines the individual EKF yaw estimates
		float _model_weights[N_MODELS_EKFGSF];
		float _gsf_yaw; 		// yaw estimate (rad)
		float _gsf_yaw_variance; 	// variance of yaw estimate (rad^2)

	}EKFGSF_YAW_t;

#endif // !EKF_EKFGSF_YAW_H
