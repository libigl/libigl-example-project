#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/cgal/convex_hull.h>

int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  const Eigen::MatrixXd V= (Eigen::MatrixXd(4,3)<<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    1.0,0.0,0.0).finished();
  Eigen::MatrixXi F;
  igl::copyleft::cgal::convex_hull(V,F);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}
