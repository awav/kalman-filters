#include <iostream>
#include "Eigen/Cholesky"
#include "ukf.h"
#include "tools.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  int weights_init = 0.5 / (n_aug_ + lambda_);
  weights_.fill(weights_init);
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  R_ << std_radr_ * std_radr_, 0, 0,
        0, std_radphi_ * std_radphi_, 0,
        0, 0, std_radrd_ * std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt_ the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt_) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

/**
 * Augment Sigma Points
 */
void UKF::AugmentSigmaPoints() {
  x_aug_.head(n_x_) = x_;
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(n_x_ + 1, n_x_ + 1) = std_a_ * std_a_;
  P_aug_(n_x_ + 2, n_x_ + 2) = std_yawdd_ * std_yawdd_;

  const auto& L =  P_aug_.llt().matrixL();

  Xsig_aug_.col(0) = x_aug_;
  auto coeff = std::sqrt(lambda_ + n_aug_);
  for (auto i = 1; i <= n_aug_; ++i) {
    Xsig_aug_.col(i) = x_aug_ + coeff * L.col(i-1);
    Xsig_aug_.col(i + n_aug_) = x_aug_ - coeff * L.col(i-1);
  }
}

/**
 * Predict Sigma Points
 */
void UKF::PredictSigmaPoints() {
  int n_aug_x = n_aug_ * 2 + 1;
  for (auto i = 0; i < n_aug_x; ++i) {
    const double p_x      = Xsig_aug_(0, i);
    const double p_y      = Xsig_aug_(1, i);
    const double v        = Xsig_aug_(2, i);
    const double yaw      = Xsig_aug_(3, i);
    const double yawd     = Xsig_aug_(4, i);
    const double nu_a     = Xsig_aug_(5, i);
    const double nu_yawdd = Xsig_aug_(6, i);

    double px_p = 0;
    double py_p = 0;

    if (std::abs(yawd) > epsilon_) {
      px_p = p_x + v/yawd * (std::sin(yaw + yawd * dt_) - std::sin(yaw));
      py_p = p_y + v/yawd * (std::cos(yaw) - std::cos(yaw + yawd * dt_));
    } else {
        px_p = p_x + v * dt_ * std::cos(yaw);
        py_p = p_y + v * dt_ * std::sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * dt_;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5 * nu_a * dt_ * dt_ * std::cos(yaw);
    py_p = py_p + 0.5 * nu_a * dt_ * dt_ * std::sin(yaw);
    v_p = v_p + nu_a * dt_;

    yaw_p = yaw_p + 0.5 * nu_yawdd * dt_ * dt_;
    yawd_p = yawd_p + nu_yawdd * dt_;

    // write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance() {
  x_ = (Xsig_aug_.array().rowwise() * weights_.array().transpose())
          .rowwise()
          .sum();

  P_.fill(0.0);
  auto diffs = Xsig_aug_.colwise() - x_;
  for (int i = 0; i < (n_aug_ * 2 + 1); ++i) {
    double phi = diffs(3, i);
    phi = std::remainder(phi + M_PI, 2 * M_PI);
    phi = phi >= 0 ? (phi - M_PI) : (phi + M_PI);
    diffs(3, i) = alpha;
    P_ += weights_(i) * diffs.col(i) * diffs.col(i).transpose();
  }
}

void UKF::PredictRadarMeasurement() {
  Zsig_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = std::cos(yaw) * v;
    double v2 = std::sin(yaw) * v;

    Zsig_(0, i) = std::sqrt(p_x*p_x + p_y*p_y);
    Zsig_(1, i) = std::atan2(p_y, p_x);
    Zsig_(2, i) = (p_x*v1 + p_y*v2) / std::sqrt(p_x*p_x + p_y*p_y);
  }

  z_pred_ = (Zsig_.array().rowwise() * weights_.array().transpose())
              .rowwise()
              .sum();

  S.fill(0.0);
  auto diffs = Zsig_.colwise() - z_pred_;
  for (int i = 0; i < 2 * n_aug + 1; i++) {
    double phi = diffs(3, i);
    phi = std::remainder(phi + M_PI, 2 * M_PI);
    phi = phi >= 0 ? (phi - M_PI) : (phi + M_PI);
    diffs(3, i) = alpha;
    S = S + weights(i) * z_diff * z_diff.transpose();
  }

}
