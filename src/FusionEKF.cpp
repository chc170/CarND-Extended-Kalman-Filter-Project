#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


#define LASER_VAR_PX 0.0225
#define LASER_VAR_PY 0.0225
#define RADAR_VAR_RHO 0.09
#define RADAR_VAR_PHI 0.0009
#define RADAR_VAR_RHO_DOT 0.09

#define SIGMA_AX 9
#define SIGMA_AY 9

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_      = MatrixXd(3, 4);

    // measurement covariance matrix - laser
    R_laser_ << LASER_VAR_PX, 0,
                0, LASER_VAR_PY;

    // measurement covariance matrix - radar
    R_radar_ << RADAR_VAR_RHO, 0, 0,
                0, RADAR_VAR_PHI, 0,
                0, 0, RADAR_VAR_RHO_DOT;

    // Kalman filter variables
    ekf_.x_ = VectorXd(4);
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.F_ = MatrixXd(4, 4);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/

    if (!is_initialized_) {
        /**
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF: " << endl;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            float rho     = measurement_pack.raw_measurements_[0];
            float phi     = measurement_pack.raw_measurements_[1];
            float rho_dot = measurement_pack.raw_measurements_[2];

            VectorXd polar(3);
            polar << rho, phi, rho_dot;
            cout << "1\n";
            ekf_.x_ = tools.PolarToCartesian(polar);
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            float px = measurement_pack.raw_measurements_[0];
            float py = measurement_pack.raw_measurements_[1];

            cout << "2\n";
            ekf_.x_ << px, py, 0, 0;
        }

        ekf_.P_ << 1, 0, 0, 0,
                   0, 1, 0, 0,
                   0, 0, 1000, 0,
                   0, 0, 0, 1000;
        cout << "3\n";
        previous_timestamp_ = measurement_pack.timestamp_;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
       * Update the state transition matrix F according to the new elapsed time.
        - Time is measured in seconds.
       * Update the process noise covariance matrix.
       * Use sigma_ax = 9 and sigma_ay = 9 for your Q matrix.
     */

    // elapsed time (dt)
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    float dt2 = dt * dt;
    float dt3 = dt * dt2;
    float dt4 = dt * dt3;

    // Update state transition matrix (F) with dt
    ekf_.F_ <<  1,  0, dt,  0,
                0,  1,  0, dt,
                0,  0,  1,  0,
                0,  0,  0,  1;

    // Updte covariance matrix (Q) with control matrix and dt
    ekf_.Q_ << dt4/4 * SIGMA_AX, 0, dt3/2 * SIGMA_AX, 0,
               0, dt4/4 * SIGMA_AY, 0, dt3/2 * SIGMA_AY,
               dt3/2 * SIGMA_AX, 0, dt2 * SIGMA_AX, 0,
               0, dt3/2 * SIGMA_AY, 0, dt2 * SIGMA_AY;

    cout << "Predict...\n";
    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        cout << "Size1: " << measurement_pack.raw_measurements_.size() << endl;
        float rho     = measurement_pack.raw_measurements_[0];
        float phi     = measurement_pack.raw_measurements_[1];
        float rho_dot = measurement_pack.raw_measurements_[2];

        VectorXd polar(3);
        polar << rho, phi, rho_dot;

        // Apply Jacobian matrix to current state
        cout << "Jacobian...\n";
        Hj_ = tools.CalculateJacobian(ekf_.x_);

        // Update Kalman filter variables
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;

        cout << "Update EKF...\n";
        // Update
        ekf_.UpdateEKF(polar);
    } else {
        // Laser updates
        cout << "Size2: " << measurement_pack.raw_measurements_.size() << endl;
        float px = measurement_pack.raw_measurements_[0];
        float py = measurement_pack.raw_measurements_[1];

        VectorXd z(2);
        z << px, py;

        // Update Kalman filter variables
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;

        cout << "Update...\n";
        // Update
        ekf_.Update(z);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
