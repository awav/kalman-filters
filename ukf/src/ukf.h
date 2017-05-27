#ifndef UKF_H
#define UKF_H

#include <vector>
#include <string>
#include <fstream>
#include <cstdint>
#include "Eigen/Dense"
#include "tools.h"
#include "measurement_package.h"

class UKF {
public:

  ///* State dimension
  const int n_x_ = 5;

  ///* Augmented state dimension
  const int n_aug_ = n_x_ + 2;

  ///* Measurement vector size
  const int n_z_ = 3;

  ///* Sigma point spreading parameter
  const double lambda_ = 3 - n_aug_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  const double std_a_ = 0.2;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  const double std_yawdd_ = 0.2;

  ///* Laser measurement noise standard deviation position1 in m
  const double std_laspx_ = 0.15;

  ///* Laser measurement noise standard deviation position2 in m
  const double std_laspy_ = 0.15;

  ///* Radar measurement noise standard deviation radius in m
  const double std_radr_ = 0.3;

  ///* Radar measurement noise standard deviation angle in rad
  const double std_radphi_ = 0.03;

  ///* Radar measurement noise standard deviation radius change in m/s
  const double std_radrd_ = 0.3;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_ = false;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_ = true;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_ = true;

  ///* Weights of sigma points
  Eigen::VectorXd weights_ = Eigen::MatrixXd(n_x_);

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_ = Eigen::VectorXd::Zero(n_x_);

  ///* state covariance matrix
  Eigen::MatrixXd P_ = Eigen::MatrixXd::Identity(n_x_, n_x_);

  ///*
  Eigen::MatrixXd R_ = Eigen::MatrixXd(n_z_, n_z_);

  ///*
  Eigen::MatrixXd S_ = Eigen::MatrixXd(n_z_, n_z_);

  ///* predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_ = Eigen::MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  ///*
  Eigen::VectorXd Zsig_ = Eigen::MatrixXd::Zero(n_z_, 2 * n_aug_ + 1);

  ///*
  Eigen::VectorXd z_pred_ = Eigen::VectorXd::Zero(n_z_);

  ///* augmented state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in
  ///* SI units and rad
  Eigen::VectorXd x_aug_ = Eigen::VectorXd::Zero(n_aug_);

  ///* augmented state covariance matrix
  Eigen::MatrixXd P_aug_ = Eigen::MatrixXd::Identity(n_aug_, n_aug_);

  ///* augmented sigma points matrix
  Eigen::MatrixXd Xsig_aug_ = Eigen::MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);

  ///* time when the state is true, in us
  int64_t time_us_ = 0;

  double dt_ = 0.0;

  const epsilon_ = 1e-7;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

 private:
  /**
   * Augment sigma points
   */
  void UKF::AugmentSigmaPoints();

  /**
   * Predict Sigma Points
   */
  void UKF::PredictSigmaPoints();
};

#endif /* UKF_H */
