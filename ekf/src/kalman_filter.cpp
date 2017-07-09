#include <iostream>
#include "kalman_filter.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;

void KalmanFilter::Init(const VectorXd &x_in, const MatrixXd &P_in) {
  x_ = x_in;
  P_ = P_in;
}


const VectorXd& KalmanFilter::StateVector() const {
  return x_;
}

const MatrixXd& KalmanFilter::CovarianceVector() const {
  return P_;
}

void KalmanFilter::Predict(const KalmanFilter::FFunc &prediction,
                           const KalmanFilter::FPrimeFunc &prediction_prime,
                           const MatrixXd &q) {
  x_ = prediction(x_);
  MatrixXd p_prime = prediction_prime();
  MatrixXd p_prime_t = p_prime.transpose();
  P_ = (p_prime * P_ * p_prime_t) + q;
}


void KalmanFilter::Update(const KalmanFilter::HFunc &h_func,
                          const KalmanFilter::HPrimeFunc &h_prime_func,
                          const MatrixXd &r) {
  VectorXd y = h_func(x_);
  MatrixXd h = h_prime_func(x_);
  MatrixXd h_t =  h.transpose();
  MatrixXd s = (h * P_ * h_t) + r;
  MatrixXd s_i = s.inverse();
  MatrixXd k = P_ * h_t * s_i;
  x_ = x_ + (k * y);
  auto l = x_.size();
  P_ = (MatrixXd::Identity(l, l) - (k * h)) * P_;
}
