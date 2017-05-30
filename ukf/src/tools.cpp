#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  auto length = estimations.size();
  VectorXd residual = VectorXd(4).setZero();

	if (length != ground_truth.size() || length == 0) {
		std::cout << "Invalid estimation or ground_truth data" << std::endl;
		return residual;
	}

  VectorXd squared_err = VectorXd(4).setZero();
  for (int i = 0; i < length; ++i) {
    residual = (estimations[i] - ground_truth[i]).array().square();
    squared_err += residual;
  }

  const VectorXd &mse = squared_err / length;
  return mse.array().sqrt();
}
