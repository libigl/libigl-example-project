// #include "remesh_botsch.h"
// #include <igl/read_triangle_mesh.h>
// #include <igl/write_triangle_mesh.h>

// using namespace std;

// int main(int argc, const char* argv[]){


//     Eigen::MatrixXd V;
//     Eigen::VectorXd target;
//     Eigen::MatrixXi F;
//     string in,out;
//     bool project = false;
//     int iterations = 10;
//     double h = 0.005;
//     in = "/home/anandhu/Desktop/test/libigl-example-project/after-out34.obj";
// 	out="/home/anandhu/Desktop/test/libigl-example-project/out.obj";
// //     if(argc>1){
// // 	    in = argv[1];
// // 	    argindex = 1;
// //     }else{
// //         std::cout<<R"(
// // No input specified. To run in command line, issue

// // ./remeshmesh [input.ext] [output.ext] [-i num_iterations] [-h target_edge_length]

// // The mesh in input.ext must be closed and manifold,
// // and ext can be either of obj, ply, off, stl or mesh.
// // )";
// //         return argindex;
// //     }
// // 	    while(argindex+1<argc){
// // 		    if(strncmp(argv[argindex+1],"-i",2)==0){
// // 			    iterations = atoi(argv[argindex+2]);
// // 		            argindex = argindex+2;
// // 		    }else{
// // 		    if(strncmp(argv[argindex+1],"-h",2)==0){
// // 			    h = atof(argv[argindex+2]);
// // 		            argindex = argindex+2;
// // 		    }else{
// // 		    if(strncmp(argv[argindex+1],"-p",2)==0){
// // 			    project = true;
// // 		            argindex = argindex+1;
// // 		    }else{
// // 		    	    out = argv[argindex+1];
// // 		            argindex = argindex+1;
// // 		    }}}
// // 	    }
//     igl::read_triangle_mesh(in,V,F);
// 	Eigen::VectorXi feature(5); // Initialize the size to 5
//     feature << 1, 2, 3, 4, 5;
//     target = Eigen::VectorXd::Constant(V.rows(),h);
//     remesh_botsch(V,F,target,iterations,feature,project);
//     igl::write_triangle_mesh(out,V,F);
// }
