#include <igl/read_triangle_mesh.h>
#include <igl/loop.h>
#include <igl/upsample.h>
#include <igl/writeOFF.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <vector>
#include "wt.h"

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace igl;
  Eigen::MatrixXi F, OF;
  Eigen::MatrixXd V, OV;

  // igl::read_triangle_mesh("/home/michelle/Documents/LIBIGL/hackathon/libigl-example-project/knightloop.off",OV,OF);
  string repo_path = "/home/michelle/Documents/LIBIGL/LABELS_BUG/mishfork/libigl-example-project/";
  repo_path += "3rdparty/libigl/tutorial/data/";
  const string mesh_path = repo_path + "decimated-knight.off";
  igl::read_triangle_mesh(mesh_path,OV,OF);

  V = OV;
  F = OF;

  // cout<<R"(Usage:
  //   1  Restore Original mesh
  //   2  Apply In-plane upsampled mesh
  //   3  Apply Loop subdivided mesh
  //   4  Apply False barycentric subdivision
  // )";
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);

  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
    {
      switch(key)
      {
        default:
          return false;
        case '1':
        {
          V = OV;
          F = OF;
          break;
        }
        case '2':
        {
          igl::upsample( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
          break;
        }
        case '3':
        {
          igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
          cout<<"Number of vertices is: " << V.rows() << endl;
          cout<<"-------"<<endl;
          // igl::writeOFF("/home/michelle/Documents/LIBIGL/hackathon/libigl-example-project/knightloop.off", V, F);
          break;
        }
        case '4':
        {
          igl::false_barycentric_subdivision(
            Eigen::MatrixXd(V),Eigen::MatrixXi(F),V,F);
          break;
        }
      }
      viewer.data().clear();
      viewer.data().set_mesh(V,F);
      viewer.data().set_face_based(true);
      return true;
    };
  // viewer.launch();


  // Begin wavelet
  std::cout << "verts: " << V.rows() << std::endl;
  std::cout << "faces: " << F.rows() << std::endl;
  
  Eigen::MatrixXi FC;
  covering_mesh(F, FC);

}
