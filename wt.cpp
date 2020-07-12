#include "wt.h"
#include <iostream>
#include <vector>
#include <map>

bool test_covering_mesh(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& V
){
  // Begin wavelet
  std::cout << "Num verts: " << V.rows() << std::endl;
  std::cout << "Num faces: " << F.rows() << std::endl;
  
  Eigen::MatrixXi F_c;
	std::vector<std::vector<int>> sub_meshes;
  covering_mesh(F, F_c);
	std::cout << "Completed covering mesh" << std::endl;
	std::cout << "Begin connected components" << std::endl;
  connected_components(F_c, sub_meshes);
	
	for(auto it=sub_meshes.begin(); it!=sub_meshes.end(); it++)
	{
		std::cout << it->size() << std::endl; 
		if(it->size()*4==F.rows())
		{
			is_equivalence(F,V,*it);
		}
	}

	return true;
};

void edge_incident_faces(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces
){
	for(int f=0; f<F.rows(); f++)
	{
		int v1 = F(f,0);
		int v2 = F(f,1);
		int v3 = F(f,2);

		// std::cout << "v1: " << v1 << std::endl;
    // std::cout << "v2: " << v2 << std::endl;
    // std::cout << "v3: " << v3 << std::endl;

		incident_faces[std::make_pair(std::min(v1,v2),std::max(v1,v2))].push_back(f);
		incident_faces[std::make_pair(std::min(v2,v3),std::max(v2,v3))].push_back(f);
		incident_faces[std::make_pair(std::min(v3,v1),std::max(v3,v1))].push_back(f);
	}
};

void covering_mesh(
	const Eigen::MatrixXi& F,
 	Eigen::MatrixXi& F_c
){
	std::map<std::pair<int,int>, std::vector<int>> incident_faces;
  edge_incident_faces(F, incident_faces);
  std::cout << "Num edges: " << incident_faces.size() << std::endl;

	std::vector<std::tuple<int, int, int> > tiles;
	for(int f=0; f<F.rows(); f++)
	{
		std::cout << "-----FACE: " << f << std::endl;
		int v1 = F(f,0);
		int v2 = F(f,1);
		int v3 = F(f,2);

		std::cout << "v1: " << v1 << std:: endl; 
		std::cout << "v2: " << v2 << std:: endl; 
		std::cout << "v3: " << v3 << std:: endl; 


		int e1v1 = std::min(v1,v2);
		int e1v2 = std::max(v1,v2);
		int e2v1 = std::min(v2,v3);
		int e2v2 = std::max(v2,v3);
		int e3v1 = std::min(v3,v1);
		int e3v2 = std::max(v3,v1);

		std::pair<int, int> e1(e1v1, e1v2);
		std::pair<int, int> e2(e2v1, e2v2);
		std::pair<int, int> e3(e3v1, e3v2);

		bool isRegular = 
			   incident_faces[e1].size()==2
			&& incident_faces[e2].size()==2
			&& incident_faces[e3].size()==2;

		if(isRegular)
		{
			std::cout << "e1: " << incident_faces[e1][0] << " " << incident_faces[e1][1] << std::endl;
			std::cout << "e2: " << incident_faces[e2][0] << " " << incident_faces[e2][1] << std::endl;
			std::cout << "e3: " << incident_faces[e3][0] << " " << incident_faces[e3][1] << std::endl;
			// Neighbouring triangles to current triangle
			int nt1 = incident_faces[e1][0]==f ? incident_faces[e1][1] : incident_faces[e1][0];
			int nt2 = incident_faces[e2][0]==f ? incident_faces[e2][1] : incident_faces[e2][0];
			int nt3 = incident_faces[e3][0]==f ? incident_faces[e3][1] : incident_faces[e3][0];

			// Get tile vertices
			int v1p, v2p, v3p;
			for(int i=0; i<3; i++)
			{
				if(F(nt1, i)!=e1v1 && F(nt1, i)!=e1v2) v1p = F(nt1, i);
				if(F(nt2, i)!=e2v1 && F(nt2, i)!=e2v2) v2p = F(nt2, i);
				if(F(nt3, i)!=e3v1 && F(nt3, i)!=e3v2) v3p = F(nt3, i);
			}
			tiles.push_back(std::tuple<int, int, int>(v1p,v2p,v3p));
			std::cout << "nt1: " << nt1 << std::endl;
			std::cout << "nt2: " << nt2 << std::endl;
			std::cout << "nt3: " << nt3 << std::endl;
		}

		F_c.setIdentity(tiles.size(),3);
		for(int t=0; t < tiles.size(); t++)
		{
			F_c(t, 0) = std::get<0>(tiles[t]);
			F_c(t, 1) = std::get<1>(tiles[t]);
			F_c(t, 2) = std::get<2>(tiles[t]);
		}
	}
	std::cout << "Num tiles: " << F_c.rows() << std::endl;
};

void connected_components(
	const Eigen::MatrixXi& tiles,
 	std::vector<std::vector<int>>& sub_meshes
){
	std::map<int, std::vector<int>*> where_are_you; // Which partition is the face in
	for(int tid=0; tid<tiles.rows(); tid++)
	{
		// At first, the partition of the mesh is 
		// a collection of singletons of the fids
		where_are_you[tid] = new std::vector<int>;
		where_are_you[tid]->emplace_back(tid);
	}
	assert(where_are_you.size()==tiles.rows());

	// Get all the edges in the mesh
	std::map<std::pair<int,int>, std::vector<int>> incident_tiles;
  edge_incident_faces(tiles, incident_tiles);

	std::map<std::pair<int,int>, std::vector<int>>::iterator it = incident_tiles.begin();
	while (it != incident_tiles.end())
	{
		if(it->second.size()==3)
		{
			std::cout << it->second[0] << std::endl;
			std::cout << it->second[1] << std::endl;
			std::cout << it->second[2] << std::endl;
			std::cout << "-------" << std::endl;
		}
		assert((it->second).size()<=2);
		if((it->second).size()==2)
		{
			int tid1 = it->second[0];
			int tid2 = it->second[1];
			if(where_are_you[tid1]!=where_are_you[tid2])
			{
				where_are_you[tid1]->insert(
					where_are_you[tid1]->end(),
					where_are_you[tid2]->begin(),
					where_are_you[tid2]->end()
				);
				delete where_are_you[tid2];
				where_are_you[tid2] = where_are_you[tid1];
			}
			else
			{
				assert(where_are_you[tid1]->size() == where_are_you[tid2]->size());
			}
		}
		it++;
	}

	std::cout << "Find unique submeshes" << std::endl;
	for(int tid=0; tid<tiles.rows(); tid++)
	{
		if(where_are_you[tid] != NULL)
		{
			std::vector<int>* temp = where_are_you[tid];
			sub_meshes.emplace_back(*temp);
			for(size_t i=0; i<temp->size(); i++)
			{
				where_are_you[ temp->at(i) ] = NULL;
			}
			delete temp;
		}
	}

	// std::cout << sub_meshes.size() << std::endl;
	// std::cout << sub_meshes.at(0).size() << std::endl;
	// std::cout << sub_meshes.at(1).size() << std::endl;
	// std::cout << sub_meshes.at(2).size() << std::endl;
	// std::cout << sub_meshes.at(3).size() << std::endl;
};

bool is_equivalence(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& V,
	const std::vector<int>& sub_mesh
){
	return true;
}

