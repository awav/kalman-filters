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
  void ProcessMeasurement(const MeasurementPackage &mpkg);

  /**
   * Getter of EKF state vector
   */
  const Eigen::VectorXd& StateVector() const;

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;

private:
  void Initialise(const MeasurementPackage &mpkg);
  void Predict(const MeasurementPackage &mpkg);
  void UpdateRadar(const MeasurementPackage &mpkg);
  void UpdateLaser(const MeasurementPackage &mpkg);

  double Eps(double val) const;

private:
  // tool object used to compute Jacobian and RMSE
  Eigen::MatrixXd F_;
  Eigen::MatrixXd Q_;
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd Hj_;

  long previous_timestamp_ = -1;
  const double epsilon_ = 1e-5;
  const double ax_ = 9;
  const double ay_ = 9;

};

#endif /* FUSION_EKF_H_ */
