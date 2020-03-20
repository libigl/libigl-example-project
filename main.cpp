#include <igl/read_triangle_mesh.h>
#include <igl/loop.h>
#include <igl/upsample.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace igl;
  Eigen::MatrixXi F;
  Eigen::MatrixXd V;
  // bool show_swept_volume = false;
  // read_triangle_mesh(
  //     TUTORIAL_SHARED_PATH "/decimated-knight.off",OV,OF);
  
  // // Inline mesh of a cube
  // const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
  //   0.0,0.0,0.0,
  //   0.0,0.0,1.0,
  //   0.0,1.0,0.0,
  //   0.0,1.0,1.0,
  //   1.0,0.0,0.0,
  //   1.0,0.0,1.0,
  //   1.0,1.0,0.0,
  //   1.0,1.0,1.0).finished();
  // const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
  //   1,7,5,
  //   1,3,7,
  //   1,4,3,
  //   1,2,4,
  //   3,8,7,
  //   3,4,8,
  //   5,7,8,
  //   5,8,6,
  //   1,5,6,
  //   1,6,2,
  //   2,6,8,
  //   2,8,4).finished().array()-1; // Test

  // Inline mesh of a tetrahedron
  const Eigen::MatrixXd OV= (Eigen::MatrixXd(4,3)<<
    0.0,0.0,0.0,
    1.0,0.0,0.0,
    0.66,0.33,1.0,
    1.0,1.0,0.0).finished();
  const Eigen::MatrixXi OF = (Eigen::MatrixXi(4,3)<<
    1,2,3,
    1,3,4,
    1,4,2,
    2,4,3).finished().array()-1; // Test

  const Eigen::MatrixXd P0 = (Eigen::MatrixXd(10,4)<<
    1.0,0.0,0.0,0.0,
    0.0,1.0,0.0,0.0,
    0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,1.0,
    .5, .5, 0., 0.,
    0.,.5,.5,0.,
    .5, 0., .5, 0,
    0.,0.,.5,.5,
    0.5, 0.0,0.,0.5,
    0., .5, 0.,.5).finished().array(); // Test 

  V = OV;
  F = OF;

  cout<<R"(Usage:
    1  Restore Original mesh
    2  Apply In-plane upsampled mesh
    3  Apply Loop subdivided mesh
    4  Apply False barycentric subdivision
  )";
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
          cout<<(P0*OV)<<endl;
          cout<<"-------"<<endl;
          V = P0*OV;
          cout<<(V)<<endl;
          break;
        }
        case '3':
        {
          igl::loop( Eigen::MatrixXd(V), Eigen::MatrixXi(F), V,F);
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
  viewer.launch();
}
