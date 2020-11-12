#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include "cubic_style_data.h"
#include "cubic_style_precomputation.h"
#include "cubic_style_single_iteration.h"
#include <iostream>
#include <string>
using namespace std;
using namespace Eigen;

MatrixXd get_cube_style(MatrixXd origV, MatrixXi F, double m_lambda){
  MatrixXd V = origV;
  cubic_style_data data;
  data.V = origV;
  data.F = F;
  data.lambda = m_lambda;
  cubic_style_precomputation(origV, F, data);
  
  int maxIters = 10; 
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
		 if(maxRelativeDiff < 0.001){ // stop condition
			 	break;
		 }
	}
  cout << "Total iters: " << iters << endl; //debug stop condition
  return V;
}



int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  
  // Read in arguments
  igl::read_triangle_mesh(argc>1 ? argv[1] : "../data/bunny.off", V, F);
  static double lambda = argc > 2 ? stod(argv[2]) : 0.1;
  
  MatrixXd cubeV = get_cube_style(V, F, lambda);
  
  igl::opengl::glfw::Viewer viewer;
  
  std::cout<<R"(
S,s      Increase, reduce lambda ('cubeness')
)";
  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    double d_lambda = 0.2;
    switch(key)
    {
      
      case 'S':
      case 's':
        {
        cout << "Test key press s\n";
        lambda = max(0.0, lambda + (key == 's' ? -d_lambda : d_lambda));
        cout << "lambda=" << lambda << endl;
        cubeV = get_cube_style(V, F, lambda);
        viewer.data().set_mesh(cubeV, F);
        }
        break;
      default:
        return false;
    }
    
    return true;
  };
  viewer.data().set_mesh(cubeV, F);
  viewer.launch();
  return 0;
}

