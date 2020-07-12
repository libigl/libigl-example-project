#ifndef WT
#define WT
#endif

#include <Eigen/Core>
#include <vector>
#include <map>

// units tests
bool test_covering_mesh(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& V
);

// WT methods
void edge_incident_faces(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces
);
void covering_mesh(
	const Eigen::MatrixXi& F,
 	Eigen::MatrixXi& F_c
);
// sub_meshes: vector of vectors containing fids 
// in a single connected component found 
void connected_components(
	const Eigen::MatrixXi& F,
 	std::vector<std::vector<int>>& sub_meshes
);
bool is_equivalence(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& V,
	const std::vector<int>& sub_mesh
);
bool is_quadrisection();
