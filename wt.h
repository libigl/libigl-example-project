#ifndef WT
#define WT

#include <Eigen/Core>
#include <vector>

void printstuff();

template <typename DerivedF>
bool edge_incident_faces(
	const Eigen::MatrixBase<DerivedF>& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces);

template <typename DerivedF>
bool covering_mesh(
	const Eigen::MatrixBase<DerivedF>& F,
  Eigen::MatrixBase<DerivedF>& FC);

#ifndef IGL_STATIC_LIBRARY
#  include "wt.cpp"
#endif
#endif

// Inline mesh of a tetrahedron
// const Eigen::MatrixXd OV= (Eigen::MatrixXd(4,3)<<
//   0.0,0.0,0.0,
//   1.0,0.0,0.0,
//   0.66,0.33,1.0,
//   1.0,1.0,0.0).finished();
// const Eigen::MatrixXi OF = (Eigen::MatrixXi(4,3)<<
//   1,2,3,
//   1,3,4,
//   1,4,2,
//   2,4,3).finished().array()-1; // Test

// const Eigen::MatrixXd P0 = (Eigen::MatrixXd(10,4)<<
//   1.0,0.0,0.0,0.0,
//   0.0,1.0,0.0,0.0,
//   0.0,0.0,1.0,0.0,
//   0.0,0.0,0.0,1.0,
//   .5, .5, 0., 0.,
//   .5, 0., .5, 0,
//   0.,.5,.5,0.,
//   0., .5, 0.,.5,
//   0.,0.,.5,.5,
//   0.5, 0.0,0.,0.5).finished().array(); // Test 