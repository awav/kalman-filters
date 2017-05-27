#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include <functional>
#include "Eigen/Dense"

class KalmanFilter {
public:

  using FFunc = std::function<Eigen::MatrixXd(const Eigen::VectorXd &)>;
  using HFunc = std::function<Eigen::MatrixXd(const Eigen::VectorXd &)>;
  using HPrimeFunc = std::function<Eigen::MatrixXd(const Eigen::VectorXd &)>;
  using FPrimeFunc = std::function<Eigen::MatrixXd()>;

  /**
   * Constructor
   */
  KalmanFilter() {};
  KalmanFilter(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in) : x_(x_in), P_(P_in) {}

  /**
   * Destructor
   */
  virtual ~KalmanFilter() {};

  /**
   * Init Initializes Kalman filter
   * @param x_in First measurement vector
   * @param P_in Init covariance matrix
   */
  void Init(const Eigen::VectorXd &x_in, const Eigen::MatrixXd &P_in);

  /**
   * State Vector Getter
   */
  const Eigen::VectorXd& StateVector() const;

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict(const FFunc &f,
               const FPrimeFunc &f_prime,
               const Eigen::MatrixXd &Q);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const HFunc &h,
              const HPrimeFunc &h_prime,
              const Eigen::VectorXd &z,
              const Eigen::MatrixXd &R);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */

private:

  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

};

#endif /* KALMAN_FILTER_H_ */
