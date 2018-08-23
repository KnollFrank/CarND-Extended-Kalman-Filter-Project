#include "kalman_filter.h"
#include "tools.h"
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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

void KalmanFilter::Predict(float delta_T) {
  /**
  TODO:
    * predict the state
  */
  //Modify the F matrix so that the time is integrated
  F_(0, 2) = delta_T;
  F_(1, 3) = delta_T;

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd K = getK(H_);
  updateEstimates(K, H_, y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  Tools tools;
  MatrixXd H_j = tools.CalculateJacobian(x_);
  // TODO: apply hint Normalizing Angles in "Tips and Tricks"
  VectorXd y = z - h(x_);
  y(1) = atan2(sin(y(1)), cos(y(1)));
  MatrixXd K = getK(H_j);
  updateEstimates(K, H_j, y);
}

MatrixXd KalmanFilter::getK(const MatrixXd& H) {
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  return K;
}

void KalmanFilter::updateEstimates(const MatrixXd& K, const MatrixXd& H, const VectorXd& y) {
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
}

VectorXd KalmanFilter::h(const VectorXd x) {
  float px = x(0);
  float py = x(1);
  float vx  = x(2);
  float vy = x(3);

  float ro = sqrt(px*px + py*py);
  if(ro < 0.00001){
    // TODO: warum +, oder ganz anders?
    px += 0.001;
    py += 0.001;
    ro = sqrt(px*px + py*py);
  }
  float theta = atan2(py, px);
  float ro_dot = (px*vx + py*vy)/ro;

  VectorXd h(3);
  h << ro, theta, ro_dot;
  return h;
}
