#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /* Calculate the RMSE here. */
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size

    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        VectorXd residual = ground_truth[i]-estimations[i];
        residual = residual.array() * residual.array();

        rmse += residual;
    }

    //calculate the mean
    rmse = rmse.array() / estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    /* Calculate a Jacobian here.*/
    MatrixXd Hj(3,4);

    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    float c1 = fmax(1.0e-8, px*px+py*py);
    float c2 = fmax(1.0e-8, sqrt(c1));
    float c3 = fmax(1.0e-8, (c1*c2));

    //compute the Jacobian matrix
    Hj <<  (px/c2), (py/c2), 0, 0,
          -(py/c1), (px/c1), 0, 0,
          py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;
}

MatrixXd Tools::PolarToCartesian(const VectorXd& x) {

    float rho = x[0];
    float phi = x[1];
    float rho_dot = x[2];

    float px = rho * cos(phi);
    float py = rho * sin(phi);

    float vx = rho_dot * cos(phi);
    float vy = rho_dot * sin(phi);

    VectorXd x_c(4);
    x_c << px, py, vx, vy;

    return x_c;
}
