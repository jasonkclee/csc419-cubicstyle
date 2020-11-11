#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "cubic_style_data.h"
// Inputs: 
//	 V  #V by 3 list of vertex positions
//	 F  #F by 3 list of triangle indices into rows of V
//   data  struct containing additional input constraints
// Outputs:
//	 data struct contains all precomputed info for cubic stylization
void cubic_style_precomputation(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, cubic_style_data& data);

void printData(cubic_style_data& data);

Eigen::MatrixXd calcD(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, std::vector<int>* faces);
