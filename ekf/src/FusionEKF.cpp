/*
 * Copyright 2017 Artem Artemev
 */

#include "FusionEKF.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  // initializing matrices

  H_laser_(0, 0) = 1;
  H_laser_(1, 1) = 1;

  F_(0, 2) = 1;
  F_(1, 3) = 1;

  // measurement covariance matrix - laser
  R_laser_ <<
        0.0225,      0,
        0,      0.0225;

  // measurement covariance matrix - radar
  R_radar_ <<
        0.09,      0,    0,
        0,    0.0009,    0,
        0,         0, 0.09;
}


/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

/**
 * Convert angle to [-PI, PI] range
 */
double pi_range(double phi) {
  phi = std::remainder(phi + M_PI, 2 * M_PI);
  phi = phi >= 0 ? (phi - M_PI) : (phi + M_PI);
  return phi;
}


/**
 * Getter of EKF state vector
 */
const VectorXd& FusionEKF::StateVector() const {
  return ekf_.StateVector();
}

/**
 * Process Measurement
 */
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  if (previous_timestamp_ == -1) {
    // first measurement

    previous_timestamp_ =
      measurement_pack.timestamp_ <= 0 ? 1LL
      : measurement_pack.timestamp_;

    double px = 0;
    double py = 0;
    MatrixXd P = MatrixXd::Identity(4, 4);
    VectorXd x = VectorXd(4).setZero();
    P(2, 2) = 10;
    P(3, 3) = 10;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_[0];
      double alpha = measurement_pack.raw_measurements_[1];

      px = rho * std::cos(alpha);
      py = rho * std::sin(alpha);

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
    }

    x(0) = std::abs(px) < epsilon_ ? 1 : px;
    x(1) = std::abs(py) < epsilon_ ? 1 : py;

    ekf_.Init(x, P);
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  const double ax = 5;
  const double ay = 5;

  F_(0, 2) = dt;
  F_(1, 3) = dt;

  Q_(0, 0) = dt_4 * ax * 0.25;
  Q_(0, 2) = dt_3 * ax * 0.5;
  Q_(1, 1) = dt_4 * ay * 0.25;
  Q_(1, 3) = dt_3 * ay * 0.5;
  Q_(2, 0) = dt_3 * ax * 0.5;
  Q_(2, 2) = dt_2 * ax;
  Q_(3, 1) = dt_3 * ay * 0.5;
  Q_(3, 3) = dt_2 * ay;

  const auto &prediction_func = [this](const auto &x) {
    return F_ * x;
  };

  const auto &prediction_prime_func = [this]() {
    return F_;
  };

  ekf_.Predict(prediction_func, prediction_prime_func, Q_);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
  /*****************************************************************************
   *  Update Radar
   ****************************************************************************/

    const auto &z = measurement_pack.raw_measurements_;

    const auto &measurement_func = [this](const auto &x) {
      double px = x[0];
      double py = x[1];
      double vx = x[2];
      double vy = x[3];

      double rho = std::sqrt(px*px + py*py);
      rho = std::abs(rho) < epsilon_ ? epsilon_ : rho;
      double alpha = pi_range(std::atan2(py, px));
      double rho_dot = (px*vx + py*vy) / rho;

      VectorXd h = VectorXd(3);
      h << rho, alpha, rho_dot;
      return h;
    };

    const auto &measurement_prime_func = [this](const auto &x) {
      double px = x[0];
      double py = x[1];
      double vx = x[2];
      double vy = x[3];

      double px_2 = px * px;
      double py_2 = py * py;

      double p2_sum = px_2 + py_2;
      double rho = std::sqrt(p2_sum);
      rho = std::abs(rho) < epsilon_ ? epsilon_ : rho;
      double p3_sum = p2_sum * rho;

      Hj_(0, 0) = px/rho;
      Hj_(0, 1) = py/rho;
      Hj_(1, 0) = -py/p2_sum;
      Hj_(1, 1) = px/p2_sum;
      Hj_(2, 0) = py*(vx*py-vy*px)/p3_sum,
      Hj_(2, 1) = px*(vy*px-vx*py)/p3_sum,
      Hj_(2, 2) = Hj_(0, 0);
      Hj_(2, 3) = Hj_(0, 1);
      return Hj_;
    };

    ekf_.Update(measurement_func, measurement_prime_func, z, R_radar_);

  } else {  // Laser updates
    const auto &z = measurement_pack.raw_measurements_;

    const auto &measurement_func = [this](const auto &x) {
      return H_laser_ * x;
    };

    const auto &measurement_prime_func = [this](const auto &) {
      return H_laser_;
    };

    ekf_.Update(measurement_func, measurement_prime_func, z, R_laser_);
  }

  // print the output
  /*
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  */
}
