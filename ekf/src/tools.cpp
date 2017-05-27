#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size

  auto length = estimations.size();
  VectorXd residual = VectorXd(4).setZero();

	if (length != ground_truth.size() || length == 0) {
		std::cout << "Invalid estimation or ground_truth data" << std::endl;
		return residual;
	}

  VectorXd squared_err = VectorXd(4).setZero();
  for (int i = 0; i < length; ++i) {
    /*
    std::cout << "Estimations: \n" << estimations[i] << std::endl;
    std::cout << "GroundTruth: \n" << ground_truth[i] << std::endl;
    */
    residual = (estimations[i] - ground_truth[i]).array().square();
    squared_err += residual;
  }

  const VectorXd &mse = squared_err / length;
  return mse.array().sqrt();
}
