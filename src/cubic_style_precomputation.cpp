#include "cubic_style_precomputation.h"
#include <vector>
#include <iostream>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/arap_rhs.h>

using namespace std;
using namespace Eigen;
typedef Eigen::SparseMatrix<double>SparseMatrixd; 
void printData(cubic_style_data& data){
	
	for(int k = 0; k < data.v_to_faces->size(); k++){
		cout << "Faces: \n";
		vector<int>* faces = data.v_to_faces->at(k);
		for(int i = 0; i < faces->size(); i++){
			cout << faces->at(i) << ", ";
		}
		cout << endl;
		
		cout << "D: \n";
		cout << data.v_to_D->at(k);
		cout << endl;
		
		cout << "W: \n";
		cout << data.v_to_W->at(k);
		cout << endl;
	}
	
}

MatrixXd calcD(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, vector<int>* faces){
	MatrixXd out = MatrixXd::Zero(3, 3 * faces->size());
	for(int k = 0; k < faces->size(); k++){ // go through each face
		int face = faces->at(k);
		for(int i = 0; i < 3; i++){ // go through each of 3 edges
			int ind1 = F(face, i);
			int ind2 = F(face, ((i+1)%3));
			out.block(0, k * 3 + i, 3, 1) = (V.row(ind2) - V.row(ind1)).transpose();
		}
	}
	return out;
}

MatrixXd calcW(const SparseMatrixd& cotan, const Eigen::MatrixXi& F, vector<int>* faces){
	MatrixXd out = MatrixXd::Zero(3 * faces->size(), 3 * faces->size()); // Diagonal matrix?? of cotan weights
	int count = 0;
	for(int k = 0; k < faces->size(); k++){ // go through each face
		int face = faces->at(k);
		for(int i = 0; i < 3; i++){ // go through each of 3 edges
			int ind1 = F(face, i);
			int ind2 = F(face, ((i+1)%3));
			out(count,count) = cotan.coeff(ind1, ind2);
			count ++;
			
		}
	}
	return out;
}

// Inputs: 
//	 V  #V by 3 list of vertex positions
//	 F  #F by 3 list of triangle indices into rows of V
//   data  struct containing additional input constraints
// Outputs:
//	 data struct contains all precomputed info for cubic stylization
void cubic_style_precomputation(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, cubic_style_data& data){
	// init lists
	cout << "Init lists\n";
	data.v_to_faces = new vector<vector<int>*>();
	for(int i = 0; i < V.rows(); i++){
		data.v_to_faces->push_back(new vector<int>());
	}
	data.v_to_D = new vector<MatrixXd>();
	data.v_to_W = new vector<MatrixXd>();
	data.local_u = MatrixXd::Zero(V.rows(), 3);  // initial values??
	data.local_z = MatrixXd::Zero(V.rows(), 3);
	data.local_u.setRandom();
	data.local_z.setRandom();
	
	data.local_p = VectorXd::Ones(V.rows()) * 0.0001; // values according to paper
	data.local_T_incr = 2.0;
	data.local_T_decr = 2.0;
	data.local_mu = 10;
	
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, data.mass);
	
	cout << "get map of vertex to list of faces\n";
	// get map of vertex to list of faces
	for(int i = 0; i < F.rows(); i++){
		for(int k = 0; k < 3; k++){
			int v = F(i,k);
			data.v_to_faces->at(v)->push_back(i);
		}
	}
	
	SparseMatrixd W;
	igl::cotmatrix(V,F,W);
	W *= -1;
	data.cot = W;
	
	// Use igl ARAP (check arap_linear_block)
	igl::arap_rhs(V,F,3,igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS, data.arap);
	
	data.solver.compute(data.cot);
	if(data.solver.info() != Success){
		cout << "data.solver.compute(data.cot) failed!\n";
	}
	
	cout << "Compute Di\n";
	// Compute Di and W
	for(int i = 0; i < data.v_to_faces->size(); i++){
		data.v_to_D->push_back(calcD(V,F,data.v_to_faces->at(i)));
		data.v_to_W->push_back(calcW(W,F,data.v_to_faces->at(i)));
	}
	
	// calc normals 
	igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, data.N);
	
}


