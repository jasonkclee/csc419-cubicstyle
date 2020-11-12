#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <Eigen/SparseCholesky>
struct cubic_style_data{
	Eigen::SparseMatrix<double> arap; // arap matrix used to calc. b (from arap_rhs)
	
	
	Eigen::SparseMatrix<double> cot; // cotan matrix
	Eigen::MatrixXd V; // Vertices
	Eigen::MatrixXi F; // Faces
	Eigen::MatrixXd N; // per vertex normals
	Eigen::MatrixXd local_u; // |V| x 3 matrix, holds u vectors for local step ADMM
	Eigen::MatrixXd local_z; // |V| x 3 matrix, holds z vectors for local step ADMM
	Eigen::VectorXd local_p; // |V| size, holds p for local step ADMM
	Eigen::SparseMatrix<double> mass; // mass matrix
	
	double lambda; // parameter - higher lambda => more cube like
	double local_mu; // parameter for ADMM in local step
	double local_T_incr; // parameter for ADMM in local step
	double local_T_decr; // parameter for ADMM in local step
	
	std::vector<std::vector<int>*>* v_to_faces; // map vertex to list of faces it bleongs to 
	std::vector<Eigen::MatrixXd>* v_to_D;  // map vertex to matrix, where each row of matrix is a rim/spoke edge vector at rest
	std::vector<Eigen::MatrixXd>* v_to_W;  // map vertex to cotangent weights of neigbor vertices
	
	
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver; // solver for global step. Precomputed with cot
	
	//std::vector<std::set<int>*>* v_to_adj;  // not needed, using igl arap function instead
};
