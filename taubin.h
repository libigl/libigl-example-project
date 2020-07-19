#ifndef WT
#define WT
#endif

#include <Eigen/Core>
#include <vector>
#include <map>

// Taubin methods
bool is_quadrisection(
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXd& V_in,
	Eigen::MatrixXi& F_old,
	Eigen::MatrixXd& V_old,
	Eigen::MatrixXi& F_new,
	Eigen::MatrixXd& V_new
);
void covering_mesh(
	const Eigen::MatrixXi& F_in,
 	Eigen::MatrixXi& tiles,
 	Eigen::MatrixXi& covered_faces
);
void connected_components(
	const Eigen::MatrixXi& tiles,
 	std::vector<std::vector<int>>& sub_meshes
);
void is_equivalence(
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXd& V_in,
	const std::vector<int>& candidate,
	const Eigen::MatrixXi& tiles,
	const Eigen::MatrixXi& covered_faces,
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