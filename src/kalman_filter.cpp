#include "kalman_filter.h"
#include <iostream>


using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

VectorXd RadarCoordsFromCartesian(const VectorXd &x_state);

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
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    VectorXd z_pred = H_ * x_;

    Update(z, z_pred);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    VectorXd z_pred = RadarCoordsFromCartesian(x_);

    Update(z, z_pred);
}


void KalmanFilter::Update(const VectorXd &z, VectorXd z_pred) {
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

// helpers

VectorXd RadarCoordsFromCartesian(const VectorXd &x_state) {
    float x, y, vx, vy;
    float rho, phi, rho_dot;

    x = float(x_state[0]);
    y = float(x_state[1]);
    vx = float(x_state[2]);
    vy = float(x_state[3]);
    // read rho, and prevent 0 value (needed for division)
    rho = fmaxf(float(sqrt(x*x + y*y)), 0.00001);
    phi = float(atan2(y, x));

    rho_dot = (x * vx + y * vy) / rho;

    VectorXd result = VectorXd(3);
    result << rho, phi, rho_dot;

    return result;
}
