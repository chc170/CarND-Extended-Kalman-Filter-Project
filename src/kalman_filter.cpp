#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

// Predict the state
void KalmanFilter::Predict() {

    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

// Update state matrix by using Kalman Filter equations
void KalmanFilter::Update(const VectorXd &z) {

    MatrixXd y = z - (H_ * x_);
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
}

// Update state matrix by using Extended Kalman Filter equations
void KalmanFilter::UpdateEKF(const VectorXd &z) {

    MatrixXd y = z - h_(x_);
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
}

// Convert Cartesian to polar coordinates
VectorXd KalmanFilter::h_(VectorXd &x) {

    VectorXd h_x(3);

    float px = x(0);
    float py = x(1);
    float vx = x(2);
    float vy = x(3);

    float rho, theta, rho_dot;

    // When the range is too small, it tends to
    // generate greater errors in bearing angle
    // which should be ignore by setting theta
    // to 0.
    rho     = fmax(1.0e-8, sqrt(px*px + py*py));
    rho_dot = rho_dot = (px*vx + py*vy) / rho;;
    theta   = 0;

    if (rho > 1.0e-3) {
        theta   = atan(py/px);
        rho_dot = (px*vx + py*vy) / rho;
    }

    h_x << rho, theta, rho_dot;
    return h_x;
}
