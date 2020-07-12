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

void sort3(int arr[]);
// WT methods
void edge_incident_faces(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces
);
void covering_mesh(
	const Eigen::MatrixXi& F,
 	Eigen::MatrixXi& F_c,
	std::map<int, std::vector<int>>& tile_sets
);
// sub_meshes: vector of vectors containing fids 
// in a single connected component found 
void connected_components(
	const Eigen::MatrixXi& F,
 	std::vector<std::vector<int>>& sub_meshes
);
void is_equivalence(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& V,
	const std::vector<int>& candidate,
	const Eigen::MatrixXi& F_c
);
bool is_quadrisection();
