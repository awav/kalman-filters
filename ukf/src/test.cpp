#include <vector>
#include <string>
#include <iostream>
#include <cstdint>
#include "Eigen/Dense"
#include "tools.h"
#include "measurement_package.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

int main() {
  MatrixXd m(2, 4);
  VectorXd v(2);
  VectorXd v2(4);
  VectorXd v3(4);

  m << 1, 23, 6, 9,
       3, 11, 7, 2;

  v << 2, 3;
  v2 << 1, 10, 100, 1000;
  v3 << 1, 10, 100, 1000;
  auto d = 1.01;

  auto out = m.colwise() - v;
  MatrixXd out2 = (m.array().rowwise() * v2.array().transpose());
  MatrixXd out3 = (m.array().rowwise() * v2.array().transpose()).rowwise().sum();
  std::cout << "+++" << std::endl;
  std::cout << out << std::endl;
  std::cout << "+++" << std::endl;
  std::cout << out2 << std::endl;
  std::cout << "+++" << std::endl;
  std::cout << out3 << std::endl;
  std::cout << "+++" << std::endl;
  std::cout << (d * v2 * v2.transpose()) << std::endl;

  /*
  MatrixXd mat = MatrixXd::Ones(2, 3);
  std::cout << mat << std::endl;
  VectorXd vec(3);
  vec << 10, 20, 30;
  std::cout << (mat.rowwise() * vec) << std::endl;
  */

  return 0;
}
