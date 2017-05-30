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
  void ProcessMeasurement(const MeasurementPackage &meas_package);

  /**
   * Getter of state vector
   */
  const Eigen::VectorXd& State() const;

private:

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
  void UpdateLaser(const MeasurementPackage &meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const MeasurementPackage &meas_package);

  /**
   * Augment sigma points
   */
  void AugmentSigmaPoints();

  /**
   * Predict Sigma Points
   */
  void PredictSigmaPoints(double dt_);

  /**
   * Predict mean and convariance
   */
  void PredictMeanAndCovariance();

  /**
   * Predict Radar Measurement
   */
  void PredictRadarMeasurement();

  /**
   * Predict Radar Measurement
   */
  void PredictLaserMeasurement();

  /**
   * Common Updates
   */
  double Update(const Eigen::VectorXd &dz,
                const Eigen::MatrixXd &Tc,
                const Eigen::MatrixXd &S);

private:

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_ = true;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_ = true;

  ///* State dimension
  const int n_x_ = 5;

  ///* Augmented state dimension
  const int n_aug_ = n_x_ + 2;

  ///* Measurement radar vector size
  const int n_z_rad_ = 3;

  ///* Measurement radar vector size
  const int n_z_las_ = 2;

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


  ///* predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  ///* augmented state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in
  ///* SI units and rad
  Eigen::VectorXd x_aug_;

  ///* augmented state covariance matrix
  Eigen::MatrixXd P_aug_;

  ///* augmented sigma points matrix
  Eigen::MatrixXd Xsig_aug_;


  ///*
  Eigen::MatrixXd Zsig_rad_;

  ///*
  Eigen::VectorXd z_pred_rad_;

  ///*
  Eigen::MatrixXd Rrad_;

  ///*
  Eigen::MatrixXd Srad_;

  ///*
  Eigen::MatrixXd Tc_rad_;


  ///*
  Eigen::MatrixXd Zsig_las_;

  ///*
  Eigen::VectorXd z_pred_las_;

  ///*
  Eigen::MatrixXd Rlas_;

  ///*
  Eigen::MatrixXd Slas_;

  ///*
  Eigen::MatrixXd Tc_las_;


  ///* time when the state is true, in us
  int64_t time_us_ = 0;

  const double epsilon_ = 1e-7;

public:

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  ///* Weights of sigma points
  Eigen::VectorXd weights_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  ///* state covariance matrix
  Eigen::MatrixXd P_;

};

#endif /* UKF_H */
