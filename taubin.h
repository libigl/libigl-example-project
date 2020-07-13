#ifndef WT
#define WT
#endif

#include <Eigen/Core>
#include <vector>
#include <map>

// Taubin methods
bool is_quadrisection(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& V,
	Eigen::MatrixXi& F_old,
	Eigen::MatrixXd& V_old,
	Eigen::MatrixXi& F_new,
	Eigen::MatrixXd& V_new
);
void covering_mesh(
	const Eigen::MatrixXi& F,
 	Eigen::MatrixXi& F_c,
	std::map<int, std::vector<int>>& tile_sets
);
void connected_components(
	const Eigen::MatrixXi& F,
	// sub_meshes: vector of vectors containing fids 
	// in a single connected component found
 	std::vector<std::vector<int>>& sub_meshes
);
void is_equivalence(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& V,
	const std::vector<int>& candidate,
	const Eigen::MatrixXi& F_c,
	Eigen::MatrixXi& F_old,
	Eigen::MatrixXd& V_old,
	Eigen::MatrixXi& F_new,
	Eigen::MatrixXd& V_new
);

// Helper methods
void edge_incident_faces(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces
);
void sort3(int arr[]);