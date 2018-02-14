#include <igl/viewer/Viewer.h>
#include <igl/readOFF.h>
#include <igl/cotmatrix.h>

#include <gmpxx.h>

#include<Eigen/SparseLU>
#include<Eigen/SparseCholesky>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;


typedef mpq_class myfloat;

typedef Eigen::Matrix<myfloat,Eigen::Dynamic,Eigen::Dynamic> MatrixXmp;
typedef Eigen::Matrix<myfloat,Eigen::Dynamic,1> VectorXmp;

int main(int argc, char *argv[])
{
  igl::readOFF("/Users/daniele/git/libigl/tutorial/shared/camelhead.off", V, F);
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  
  Eigen::SparseMatrix<double> id_m(L.rows(), L.cols());
  id_m.setIdentity();

  Eigen::SparseMatrix<myfloat> Lmp = id_m.cast<myfloat>() - L.cast<myfloat>();

  VectorXmp bmp = VectorXmp::Constant(Lmp.rows(),0);

  Eigen::SparseLU<Eigen::SparseMatrix<myfloat> > solver;
  solver.compute(Lmp);
  if(solver.info()!=Eigen::Success) {
    std::cerr << "Failed decomposition" << std::endl;
    exit(-1);
  }
  VectorXmp x = solver.solve(bmp);
  if(solver.info()!=Eigen::Success) {
    std::cerr << "Failed solve" << std::endl;
    exit(-1);
  }

  Eigen::VectorXd xout(x.rows());
  for(unsigned i=0;i<xout.rows();++i)
    xout(i) = x(i).get_d();
  
  std::cerr << xout << std::endl;
}
