#include <vector>
#include <iostream>
#include "Eigen/Dense"
#include <math.h>

int main() {
	Eigen::MatrixXd err(4,4);
	err << 2,2,2,2,
	       0,0,0,0,
				 2,2,2,2,
	       0,0,0,0;
  err.array().square();
	std::cout << err << std::endl;
	std::cout << "=====\n";
  const Eigen::VectorXd &mse = err.rowwise().sum() * 1/2;
	std::cout << mse << std::endl;
	std::cout << "=====\n";
	const Eigen::VectorXd &smse = mse.array().sqrt();
	std::cout << smse << std::endl;
	std::cout << "=====\n";


  std::vector<double> alphas {
		M_PI,
		M_PI - M_PI / 2,
		M_PI*2,
		M_PI*2 - M_PI / 2,
		M_PI*3,
		M_PI*3 - M_PI / 2,
		M_PI*4,
		M_PI*4 - M_PI / 2,
		M_PI*5,
		M_PI*5 - M_PI / 2,
		M_PI*6,
		M_PI*6 - M_PI / 2
	};
	for (double alpha: alphas) {
		double tmp = alpha;
    alpha = std::remainder(alpha + M_PI, 2 * M_PI);
    //alpha = alpha >= 0 ? (alpha - M_PI) : (alpha + M_PI);
		std::cout << "Before: " << tmp << ", After: " << alpha << std::endl;
	}

  return 0;
}
