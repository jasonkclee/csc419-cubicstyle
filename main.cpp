#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include "cubic_style_data.h"
#include "cubic_style_precomputation.h"
#include "cubic_style_single_iteration.h"
#include <iostream>
#include <string>
using namespace std;
using namespace Eigen;
void getTestData(Eigen::MatrixXd& V, Eigen::MatrixXi& F){
  V = MatrixXd(4,3);
  F = MatrixXi(4,3);
	V << 0,0,0, 0.5, 0, 1,  1, 0, 0,  0.5, 1, 0.5;
	F << 0, 2, 1,   0, 3, 2,  2, 3, 1,  1,3,0;
}
int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // Load in a mesh
  igl::read_triangle_mesh(argc>1 ? argv[1] : "../data/bunny.off", V, F);
  
  cubic_style_data data;
  data.lambda = argc > 2 ? stod(argv[2]) : 0.1;
  cout << "lambda: " << data.lambda;
  //getTestData(V,F);
  
  Eigen::MatrixXd origV = V;
  
  data.V = V;
  data.F = F;
  cubic_style_precomputation(V, F, data);

	int maxIters = 50; // ??
	double maxRelativeDiff; // largest change in a single iteration (value is relative to difference from undeformed vertices)
	int iters = 0;
	for(; iters < maxIters; iters++){
		 MatrixXd prevV = V;
		 cubic_style_single_iteration(data, V);
		 MatrixXd iterDiff = V - prevV;
		 MatrixXd totalDiff = V - origV;
		 double maxIterDiff = iterDiff.array().abs().maxCoeff();
		 double maxTotalDiff = totalDiff.array().abs().maxCoeff();
		 maxRelativeDiff = maxIterDiff / maxTotalDiff;
		 cout << "maxRelativeDiff: " << maxRelativeDiff << endl;
		 if(maxRelativeDiff < 0.001){ // stop condition
			 	break;
		 }
	}
  cout << "Total iters: " << iters << endl; //debug stop condition
  MatrixXd diff = V - origV;
  double norms = 0;
  for(int i =0; i < V.rows(); i++){
  	norms += diff.row(i).norm();
  }
  cout << "Average diff: " << norms / (double) V.rows() << endl;
  
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.launch();
  return 0;
}

