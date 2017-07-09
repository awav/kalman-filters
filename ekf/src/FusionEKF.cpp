/*
 * Copyright 2017 Artem Artemev
 */

#include <iostream>

#include "FusionEKF.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
    : F_(MatrixXd::Identity(4, 4)), Q_(MatrixXd::Zero(4, 4)),
      R_laser_(MatrixXd(2, 2)), R_radar_(MatrixXd(3, 3)),
      H_laser_(MatrixXd::Zero(2, 4)), Hj_(MatrixXd::Zero(3, 4)), ax_(9),
      ay_(9) {
  // initializing matrices

  H_laser_(0, 0) = 1;
  H_laser_(1, 1) = 1;

  F_(0, 2) = 1;
  F_(1, 3) = 1;

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0, 0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

/**
 * Getter of EKF state vector
 */
const VectorXd &FusionEKF::StateVector() const { return ekf_.StateVector(); }

double FusionEKF::Eps(double val) const {
  return std::abs(val) > epsilon_ ? val : epsilon_;
}

void FusionEKF::Initialise(const MeasurementPackage &mpkg) {
  previous_timestamp_ = mpkg.timestamp_ <= 0 ? 1LL : mpkg.timestamp_;
  double px = 0;
  double py = 0;
  double vx = 0;
  double vy = 0;
  MatrixXd P = MatrixXd::Identity(4, 4);
  VectorXd x = VectorXd(4).setZero();
  P(2, 2) = 1000;
  P(3, 3) = 1000;
  if (mpkg.sensor_type_ == MeasurementPackage::RADAR) {
    const double rho = mpkg.raw_measurements_[0];
    const double alpha = mpkg.raw_measurements_[1];
    const double rho_dot = mpkg.raw_measurements_[2];
    const double cos_alpha = std::cos(alpha);
    const double sin_alpha = std::sin(alpha);
    px = rho * cos_alpha;
    py = rho * sin_alpha;
    vx = rho_dot * cos_alpha;
    vy = rho_dot * sin_alpha;
  } else if (mpkg.sensor_type_ == MeasurementPackage::LASER) {
    px = mpkg.raw_measurements_[0];
    py = mpkg.raw_measurements_[1];
  }
  x(0) = std::abs(px) < epsilon_ ? 1 : px;
  x(1) = std::abs(py) < epsilon_ ? 1 : py;
  x(2) = vx;
  x(3) = vy;
  ekf_.Init(x, P);
}

void FusionEKF::Predict(const MeasurementPackage &mpkg) {
  double dt = (mpkg.timestamp_ - previous_timestamp_) / 1.0e6;
  previous_timestamp_ = mpkg.timestamp_;

  // Avoid to often updates
  if (dt < 0.001) {
    return;
  }

  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  F_(0, 2) = dt;
  F_(1, 3) = dt;

  Q_(0, 0) = dt_4 * ax_ * 0.25;
  Q_(0, 2) = dt_3 * ax_ * 0.5;
  Q_(1, 1) = dt_4 * ay_ * 0.25;
  Q_(1, 3) = dt_3 * ay_ * 0.5;
  Q_(2, 0) = dt_3 * ax_ * 0.5;
  Q_(2, 2) = dt_2 * ax_;
  Q_(3, 1) = dt_3 * ay_ * 0.5;
  Q_(3, 3) = dt_2 * ay_;

  const auto &prediction_func = [this](const auto &x) { return F_ * x; };
  const auto &prediction_prime_func = [this]() { return F_; };
  ekf_.Predict(prediction_func, prediction_prime_func, Q_);
}

void FusionEKF::UpdateRadar(const MeasurementPackage &mpkg) {
  const auto &measurement_func = [&mpkg, this](const auto &x) {
    const double px = x[0];
    const double py = x[1];
    const double vx = x[2];
    const double vy = x[3];

    const double rho = this->Eps(std::sqrt(px * px + py * py));
    const double alpha = std::atan2(py, px);
    const double rho_dot = (px * vx + py * vy) / rho;

    VectorXd h = VectorXd(3);
    h << rho, alpha, rho_dot;

    VectorXd z = mpkg.raw_measurements_;
    VectorXd y = z - h;
    y(1) = pi_range(y(1));
    return y;
  };

  const auto &measurement_prime_func = [this](const auto &x) {
    const double px = x[0];
    const double py = x[1];
    const double vx = x[2];
    const double vy = x[3];

    const double px_2 = px * px;
    const double py_2 = py * py;

    const double p2_sum = this->Eps(px_2 + py_2);
    const double rho = this->Eps(std::sqrt(p2_sum));
    const double p3_sum = this->Eps(p2_sum * rho);

    // Calculate Jacobian
    Hj_(0, 0) = px / rho;
    Hj_(0, 1) = py / rho;
    Hj_(1, 0) = -py / p2_sum;
    Hj_(1, 1) = px / p2_sum;
    Hj_(2, 0) = py * (vx * py - vy * px) / p3_sum;
    Hj_(2, 1) = px * (vy * px - vx * py) / p3_sum;
    Hj_(2, 2) = px / rho;
    Hj_(2, 3) = py / rho;
    return Hj_;
  };

  ekf_.Update(measurement_func, measurement_prime_func, R_radar_);
}


void FusionEKF::UpdateLaser(const MeasurementPackage &mpkg) {
  const auto &measurement_func = [&mpkg, this](const auto &x) {
    MatrixXd y = mpkg.raw_measurements_ - H_laser_ * x;
    return y;
  };

  const auto &measurement_prime_func = [this](const auto &) {
    return H_laser_;
  };

  ekf_.Update(measurement_func, measurement_prime_func, R_laser_);
}

  /**
   * Process Measurement
   */
void FusionEKF::ProcessMeasurement(const MeasurementPackage &mpkg) {
  if (previous_timestamp_ < 0) {
    // first measurement
    Initialise(mpkg);
    return;
  }

  if (mpkg.sensor_type_ == MeasurementPackage::LASER) {
    Predict(mpkg);
    UpdateLaser(mpkg);
  } else if (mpkg.sensor_type_ == MeasurementPackage::RADAR) {
    Predict(mpkg);
    UpdateRadar(mpkg);
  }
}
