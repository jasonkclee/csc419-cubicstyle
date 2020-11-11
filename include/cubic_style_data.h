#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <Eigen/SparseCholesky>
struct cubic_style_data{
	Eigen::SparseMatrix<double> cot; // cotan matrix
	Eigen::MatrixXd V; // Vertices
	Eigen::MatrixXi F; // Faces
	Eigen::MatrixXd N; // per vertex normals
	Eigen::MatrixXd local_u; // |V| x 3 matrix, holds u vectors for local step ADMM
	Eigen::MatrixXd local_z; // |V| x 3 matrix, holds z vectors for local step ADMM
	Eigen::VectorXd local_p; // |V| size, holds p for local step ADMM
	Eigen::VectorXd local_a; // |V| size, holds bary area for local step
	
	double lambda;
	double local_mu; // parameter for ADMM in local step
	double local_T_incr;
	double local_T_decr;
	
	std::vector<std::vector<int>*>* v_to_faces;
	std::vector<std::set<int>*>* v_to_adj;
	std::vector<Eigen::MatrixXd>* v_to_D;  
	std::vector<Eigen::MatrixXd>* v_to_W;  // not sure
	
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver; // solver for global step using cot matrix
	
	
	
	// Declare the members that can be precomputed
	// Write some brief documentation for each component
};
