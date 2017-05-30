#include <iostream>
#include "Eigen/Cholesky"
#include "ukf.h"
#include "tools.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

double pi_range(double phi) {
  phi = std::remainder(phi + M_PI, 2 * M_PI);
  phi = phi >= 0 ? (phi - M_PI) : (phi + M_PI);
  return phi;
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  const auto n = 2 * n_aug_ + 1;

  Xsig_pred_ = MatrixXd::Zero(n_x_, n);
  x_aug_     = VectorXd::Zero(n_aug_);
  P_aug_     = MatrixXd::Identity(n_aug_, n_aug_);
  Xsig_aug_  = MatrixXd::Zero(n_aug_, n);

  Zsig_rad_   = MatrixXd::Zero(n_z_rad_, n);
  z_pred_rad_ = VectorXd::Zero(n_z_rad_);
  Rrad_       = MatrixXd(n_z_rad_, n_z_rad_);
  Srad_       = MatrixXd(n_z_rad_, n_z_rad_);
  Tc_rad_     = MatrixXd(n_x_, n_z_rad_);

  Zsig_las_   = MatrixXd::Zero(n_z_las_, n);
  z_pred_las_ = VectorXd::Zero(n_z_las_);
  Rlas_       = MatrixXd(n_z_las_, n_z_las_);
  Slas_       = MatrixXd(n_z_las_, n_z_las_);
  Tc_las_     = MatrixXd(n_x_, n_z_las_);

  weights_ = VectorXd(n);
  x_       = VectorXd::Zero(n_x_);
  P_       = MatrixXd::Identity(n_x_, n_x_);

  double weights_init = 0.5 / (n_aug_ + lambda_);
  weights_.fill(weights_init);
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  auto std_radr_2 = std_radr_ * std_radr_;
  auto std_radphi_2 = std_radphi_ * std_radphi_;
  auto std_radrd_2 = std_radrd_ * std_radrd_;
  Rrad_ << std_radr_2, 0,            0,
           0,          std_radphi_2, 0,
           0,          0,            std_radrd_2;

  auto std_laspx_2 = std_laspx_*std_laspx_;
  auto std_laspy_2 = std_laspy_*std_laspy_;
  Rlas_ << std_laspx_2, 0,
           0, std_laspy_2;
}


UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage &mp) {
  if (time_us_ == 0) {
    if (mp.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = mp.raw_measurements_(0);
      double phi = mp.raw_measurements_(1);
      double rho_dot = mp.raw_measurements_(2);

      double px = rho * std::cos(phi);
      double py = rho * std::sin(phi);
      double vx = rho_dot * std::cos(phi);
      double vy = rho_dot * std::sin(phi);

      x_ << px, py, std::sqrt(vx*vx + vy*vy), 0, 0;

    } else {

      double px = mp.raw_measurements_(0);
      double py = mp.raw_measurements_(1);

      x_ << px, py, 0, 0, 0;

    }
    time_us_ = mp.timestamp_;
    return;
  }

  double dt = (mp.timestamp_ - time_us_) / 1000000.0;
  time_us_ = mp.timestamp_;

  Prediction(dt);
  if (mp.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    PredictRadarMeasurement();
    UpdateRadar(mp);
  } else if (mp.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    PredictLaserMeasurement();
    UpdateLaser(mp);
  }
}

/**
 * Getter of state
 */
const VectorXd& UKF::State() const {
  return x_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt_ the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt_) {
  AugmentSigmaPoints();
  PredictSigmaPoints(dt_);
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLaser(const MeasurementPackage &mp) {
  VectorXd dz = VectorXd(n_z_las_);
  dz << mp.raw_measurements_(0), mp.raw_measurements_(1);
  dz -= z_pred_las_;
  NIS_laser_ = Update(dz, Tc_las_, Slas_);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage &mp) {
  VectorXd dz = VectorXd(n_z_rad_);
  dz << mp.raw_measurements_(0), mp.raw_measurements_(1), mp.raw_measurements_(2);
  dz -= z_pred_rad_;
  NIS_radar_ = Update(dz, Tc_rad_, Srad_);
}

/**
 * Augment Sigma Points
 */
void UKF::AugmentSigmaPoints() {
  x_aug_ << x_.array(), 0, 0;
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(n_x_, n_x_) = std_a_ * std_a_;
  P_aug_(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

  const MatrixXd &L =  P_aug_.llt().matrixL();

  Xsig_aug_.col(0) = x_aug_;
  double coeff = std::sqrt(lambda_ + n_aug_);
  for (auto i = 1; i <= n_aug_; ++i) {
    Xsig_aug_.col(i) = x_aug_ + coeff * L.col(i-1);
    Xsig_aug_.col(i + n_aug_) = x_aug_ - coeff * L.col(i-1);
  }
}

/**
 * Predict Sigma Points
 */
void UKF::PredictSigmaPoints(double dt_) {
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
  P_.fill(0);
  x_.fill(0);

  const int n = n_aug_ * 2 + 1;

  for (int i = 0; i < n; ++i) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  MatrixXd x_diffs = Xsig_pred_.colwise() - x_;
  for (int i = 0; i < n; ++i) {
    x_diffs(3, i) = pi_range(x_diffs(3, i));
    P_ += weights_(i) * x_diffs.col(i) * x_diffs.col(i).transpose();
  }
}


void UKF::PredictRadarMeasurement() {
  Zsig_rad_.fill(0);
  Srad_.fill(0);
  Tc_rad_.fill(0);
  z_pred_rad_.fill(0);

  const int n = 2 * n_aug_ + 1;

  for (int i = 0; i < n; i++) {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    p_x = std::abs(p_x) < epsilon_ ? epsilon_ : p_x;
    p_y = std::abs(p_y) < epsilon_ ? epsilon_ : p_y;

    double v1 = std::cos(yaw) * v;
    double v2 = std::sin(yaw) * v;

    Zsig_rad_(0, i) = std::sqrt(p_x*p_x + p_y*p_y);
    Zsig_rad_(1, i) = std::atan2(p_y, p_x);
    Zsig_rad_(2, i) = (p_x*v1 + p_y*v2) / std::sqrt(p_x*p_x + p_y*p_y);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred_rad_ += weights_(i) * Zsig_rad_.col(i);
  }

  MatrixXd z_diffs = Zsig_rad_.colwise() - z_pred_rad_;
  MatrixXd x_diffs = Xsig_pred_.colwise() - x_;

  Srad_ += Rrad_;
  for (int i = 0; i < n; i++) {
    z_diffs(1, i) = pi_range(z_diffs(1, i));
    Srad_ += weights_(i) * z_diffs.col(i) * z_diffs.col(i).transpose();

    x_diffs(3, i) = pi_range(x_diffs(3, i));
    Tc_rad_ += weights_(i) * x_diffs.col(i) * z_diffs.col(i).transpose();
  }
}


void UKF::PredictLaserMeasurement() {
  Zsig_las_.fill(0);
  Slas_.fill(0);
  Tc_las_.fill(0);
  z_pred_las_.fill(0);

  const int n = 2 * n_aug_ + 1;

  for (int i = 0; i < n; i++) {
    Zsig_las_(0, i) = Xsig_pred_(0, i);
    Zsig_las_(1, i) = Xsig_pred_(1, i);
  }

  for (int i = 0; i < n; ++i) {
    z_pred_las_ += weights_(i) * Zsig_las_.col(i);
  }

  const MatrixXd &z_diffs = Zsig_las_.colwise() - z_pred_las_;
  const MatrixXd &x_diffs = Xsig_pred_.colwise() - x_;

  Slas_ += Rlas_;
  for (int i = 0; i < n; i++) {
    Slas_ += (weights_(i) * z_diffs.col(i) * z_diffs.col(i).transpose());
    Tc_las_ += weights_(i) * x_diffs.col(i) * z_diffs.col(i).transpose();
  }
}


double UKF::Update(const VectorXd &dz, const MatrixXd &Tc, const MatrixXd &S) {
  const MatrixXd &Sinv = S.inverse();
  MatrixXd K = Tc * Sinv;

  x_ = x_ + (K * dz);
  P_ = P_ - (K * S * K.transpose());
  double nis = dz.transpose() * Sinv * dz;
  return nis;
}
