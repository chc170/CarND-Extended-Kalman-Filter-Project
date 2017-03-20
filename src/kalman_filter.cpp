#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

void KalmanFilter::Predict() {
    /* predict the state */
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /* update the state by using Kalman Filter equations */
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    MatrixXd y = z - H_ * x_;
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();

    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /* update the state by using Extended Kalman Filter equations */
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    MatrixXd y = z - h(x_);
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();

    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
}

MatrixXd KalmanFilter::h(MatrixXd &x) {

	MatrixXd h_x(3);

    float px = x(0);
	float py = x(1);
	float vx = x(2);
	float vy = x(3);

    float rho, theta, rho_dot;

    rho = sqrt(px*px + py*py);
    if (fabs(rho) < 0.0001) {
        count << "Division by zero" << endl;
        return h_x;
    }
    theta = atan(py/px);
    rho_dot = (px*vx + py*vy) / rho;
    h_x << rho, theta, rho_dot;

    return h_x;

}
