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
                        MatrixXd &H_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict(float dt, float noise_ax, float noise_ay) {
  /**
  TODO:
    * predict the state
  */
  updateF(dt);
  updateQ(dt, noise_ax, noise_ay);

  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::updateF(float dt) {
  //Modify the F matrix so that the time is integrated
  F_(0, 2) = dt;
  F_(1, 3) = dt;
}

void KalmanFilter::updateQ(float dt, float noise_ax, float noise_ay) {
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  //set the process covariance matrix Q
  Q_ << dt_4 / 4 * noise_ax, 0,                   dt_3 / 2 * noise_ax, 0,
        0,                   dt_4 / 4 * noise_ay, 0,                   dt_3 / 2 * noise_ay,
        dt_3 / 2 * noise_ax, 0,                   dt_2 * noise_ax,     0,
        0,                   dt_3 / 2 * noise_ay, 0,                   dt_2 * noise_ay;
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd &R) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  updateEstimates(getK(H_, R), H_, y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd &R) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd y = z - h(x_);
  y(1) = normalizeAngle(y(1));

  Tools tools;
  MatrixXd H_j = tools.CalculateJacobian(x_);

  updateEstimates(getK(H_j, R), H_j, y);
}

float KalmanFilter::normalizeAngle(float angle)
{
  return atan2(sin(angle), cos(angle));
}

MatrixXd KalmanFilter::getK(const MatrixXd &H, const MatrixXd &R) {
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd K = P_ * Ht * S.inverse();
  return K;
}

void KalmanFilter::updateEstimates(const MatrixXd &K, const MatrixXd &H, const VectorXd &y) {
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
}

VectorXd KalmanFilter::h(const VectorXd &x) {
  float px = x(0);
  float py = x(1);
  float vx = x(2);
  float vy = x(3);

  float ro = sqrt(px*px + py*py);
  float theta = atan2(py, px);
  float ro_dot = (px*vx + py*vy)/ro;

  VectorXd h(3);
  h << ro, theta, ro_dot;
  return h;
}
