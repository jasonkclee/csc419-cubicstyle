#include "cubic_style_single_iteration.h"
#include <iostream>
#include "cubic_style_precomputation.h"
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <algorithm>

using namespace Eigen;
using namespace std;
void print_dim(MatrixXd m){
	cout << m.rows() << ", " << m.cols() << endl;
}

// Calc rotation for a single vertex
Matrix3d calc_rot(MatrixXd D, Vector3d n, MatrixXd W, double p, MatrixXd defD, Vector3d z, Vector3d u){
	int dim = D.cols() + 1;
	MatrixXd M1(3, dim);
	M1 << D, n;
	MatrixXd M2 = MatrixXd::Zero(dim, dim);
	M2.block(0,0, dim-1, dim-1) = W;
	M2(dim-1, dim-1) = p;
	MatrixXd M3 = MatrixXd::Zero(dim, 3);
	Vector3d diff = z - u;
	M3 << defD.transpose(), diff.transpose();
	Matrix3d M = M1 * M2 * M3;
	
	JacobiSVD<Matrix3d> svd;
	svd.compute(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Matrix3d U = svd.matrixU();
	Matrix3d V = svd.matrixV();
	Matrix3d out = V * U.transpose();
	return V * U.transpose();
}

Vector3d calc_z(Matrix3d R, Vector3d n, Vector3d u, double lambda, double a, double p){
	Vector3d x = R * n + u;
	Vector3d k = Vector3d::Ones() * lambda * a / p;
	
	Vector3d out = Vector3d::Zero();
	Vector3d d1 = x - k; // Do "shrinkage" from section 4.4.3 from ADMM paper
	Vector3d d2 = -x - k;
	for(int i = 0; i < 3; i++){
		d1(i) = max(0.0, d1(i));
		d2(i) = max(0.0, d2(i));
	}
	return d1 - d2;
}

// Return set of rotations to use 
void localStep(cubic_style_data & data, Eigen::MatrixXd & U, Eigen::MatrixXd& rots){
	int maxIters = 10;
	
	
	for(int iter = 0; iter < maxIters; iter ++){
		// go through each vertex
		for(int i = 0; i < U.rows(); i++){
			
			// retrieve matrix of deformed rim/spoke edge vectors
			MatrixXd D = data.v_to_D->at(i);
			MatrixXd W = data.v_to_W->at(i);
			MatrixXd defD = calcD(U, data.F, data.v_to_faces->at(i));
			Vector3d n = data.N.row(i);
			Vector3d u = data.local_u.row(i);
			Vector3d z = data.local_z.row(i);
			double a = data.mass.coeff(i,i);
			double p = data.local_p(i);
			
			
			// calc rotation using Di, ~Di, Wi, p, ^ni , z
			Matrix3d n_R = calc_rot(D, n, W, p, defD, z, u);
			
			// calc z using n_R, ^ni, u, lambda, ai, p
			Vector3d n_z = calc_z(n_R, n, u, data.lambda, a, p);
			
			// calc u
			Vector3d n_u = u + n_R * n - n_z;
			
			// calc p
			double n_p = p;
			double r = (n_z-n_R*n).norm();
			double s = (-p * (n_z - z)).norm();
			
			// update according to 3.4.1 section referred to from ADMM paper referenced
			if(r > data.local_mu * s){
				n_p = p * data.local_T_incr;
				n_u /= data.local_T_incr; // adjust u
			}
			else if(s > data.local_mu * r){
				n_p = p / data.local_T_decr;
				n_u *= data.local_T_decr;
			}
			
			// update
			data.local_u.row(i) = n_u.transpose();
			data.local_p(i) = n_p;
			data.local_z.row(i) = n_z.transpose();
			
			if(iter == maxIters-1){
				rots.block(3*i, 0, 3,3) = n_R;
			}
		}
	}
}

// Return updated vertices
void globalStep(cubic_style_data & data, Eigen::MatrixXd & U, Eigen::MatrixXd& rots){
    // L * defU = b
	// compute b
	MatrixXd b = MatrixXd::Zero(U.rows(), 3);
	VectorXd rotCol = VectorXd::Zero(9 * U.rows());
	MatrixXd rotsT = rots;
	
	// need to flatten rots according to arap_rhs.h
	int count = 0;
	for(int i = 0; i < 3; i ++){
		for(int j = 0; j < 3; j++){
			for(int k = 0; k < U.rows(); k++){ // go through each matrix
				rotCol(count) = rots(k*3 + j, i);
				count ++;
			}
		}
	}
	
	// Check arap_rhs for rot
	VectorXd Bcol = data.arap * rotCol;
	
	/*  // NOT NEEDED, use igl arap functions instead
	for(int i = 0; i < U.rows(); i++){
		set<int>* adj = data.v_to_adj->at(i);
		Matrix3d Ri = rots.block(i*3, 0, 3, 3);
		
		for(auto it = adj->begin(); it != adj->end(); ++it){
			int v_adj = *it;
			Matrix3d Rj = rots.block(v_adj*3, 0, 3, 3);
			b.row(i) += (data.cot.coeff(i, v_adj) * 0.5 * (Ri + Rj) * (data.V.row(i) - data.V.row(v_adj)).transpose()).transpose();
			
		}
	}*/
	
	
	// solve each dim
	for(int i = 0; i < 3; i++){
		VectorXd bi = Bcol.block(i * U.rows(), 0, U.rows(), 1);
		VectorXd x = data.solver.solve(bi);
		
		if(data.solver.info() != Success){
			cout << "solver.solve(bi) failed!\n";
		}
		U.col(i) = x;
	}
	
}

// Given precomputed data and current positions of all vertices 'U', conduct a single
// iteration of the local-global solver for minimizing the cubic stylization energy.
// Output the positions of all vertices of the stylized mesh by overwriting 'U'
//
// Inputs:
// data struct that contains all the precomputed information for cubic stylization
// U #V by 3 list of current stylized mesh vertex positions
// Outputs:
// U #V by 3 list of new stylized mesh vertex positions
void cubic_style_single_iteration( cubic_style_data & data, Eigen::MatrixXd & U){
	
	MatrixXd rots(3 * U.rows(), 3);
	// do local step for rotations
	localStep(data, U, rots);
	// do global step(R) for vertices
	globalStep(data, U, rots);	
		
}
