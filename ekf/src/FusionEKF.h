/*
 * Copyright 2017 Artem Artemev
 */

#ifndef FUSION_EKF_H_
#define FUSION_EKF_H_

#include <vector>
#include <string>
#include <fstream>

#include "Eigen/Dense"
#include "measurement_package.h"
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
 public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
   * Getter of EKF state vector
   */
  const Eigen::VectorXd& StateVector() const;

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;

 private:
  // previous timestamp
  long previous_timestamp_ = -1;

  const double epsilon_ = 1e-7;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  Eigen::MatrixXd F_ = Eigen::MatrixXd::Identity(4, 4);
  Eigen::MatrixXd Q_ = Eigen::MatrixXd::Zero(4, 4);
  Eigen::MatrixXd R_laser_ = Eigen::MatrixXd(2, 2);
  Eigen::MatrixXd R_radar_ = Eigen::MatrixXd(3, 3);
  Eigen::MatrixXd H_laser_ = Eigen::MatrixXd::Zero(2, 4);
  Eigen::MatrixXd Hj_ = Eigen::MatrixXd::Zero(3, 4);
};

#endif /* FUSION_EKF_H_ */
