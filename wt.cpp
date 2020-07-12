#include "wt.h"
#include <igl/upsample.h>
#include <iostream>
#include <vector>
#include <map>

bool test_covering_mesh(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& V
){
  // Begin wavelet
  std::cout << "Num verts in input mesh: " << V.rows() << std::endl;
  std::cout << "Num faces in input mesh: " << F.rows() << std::endl;
  
  Eigen::MatrixXi F_c;
	std::vector<std::vector<int>> sub_meshes;
	std::map<int, std::vector<int>> tile_sets;
	
	std::cout << "Begin covering mesh" << std::endl;
  covering_mesh(F, F_c, tile_sets);
	std::cout << "Completed covering mesh" << std::endl;
	std::cout << "Begin connected components" << std::endl;
  connected_components(F_c, sub_meshes);
	std::cout << "Completed connected components" << std::endl;

	// Find a candidate connected component
	for(auto it=sub_meshes.begin(); it!=sub_meshes.end(); it++)
	{
		is_equivalence( F, V, *it, F_c );
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
 	Eigen::MatrixXi& F_c,
	std::map<int, std::vector<int>>& tile_sets
){
	std::map<std::pair<int,int>, std::vector<int>> incident_faces;
  edge_incident_faces(F, incident_faces);
  std::cout << "Num edges: " << incident_faces.size() << std::endl;

	std::vector<std::tuple<int, int, int> > tiles;
	for(int f=0; f<F.rows(); f++)
	{
		// std::cout << "-----FACE: " << f << std::endl;
		int v1 = F(f,0);
		int v2 = F(f,1);
		int v3 = F(f,2);

		// std::cout << "v1: " << v1 << std:: endl; 
		// std::cout << "v2: " << v2 << std:: endl; 
		// std::cout << "v3: " << v3 << std:: endl; 

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
			// std::cout << "e1: " << incident_faces[e1][0] << " " << incident_faces[e1][1] << std::endl;
			// std::cout << "e2: " << incident_faces[e2][0] << " " << incident_faces[e2][1] << std::endl;
			// std::cout << "e3: " << incident_faces[e3][0] << " " << incident_faces[e3][1] << std::endl;
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
			// std::cout << "nt1: " << nt1 << std::endl;
			// std::cout << "nt2: " << nt2 << std::endl;
			// std::cout << "nt3: " << nt3 << std::endl;
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
};

void is_equivalence(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& V,
	const std::vector<int>& candidate,
	const Eigen::MatrixXi& F_c
){
	// First test
	if(candidate.size()*4==F.rows())
	{
		std::cout << "--Begin analyzing the candidate connected component--" << std::endl;
		std::vector<int> V_i; // All the vids from og V in the tile set
		Eigen::MatrixXi submesh; // Faces in the tile set
		Eigen::MatrixXd submesh_vertices; // Vertex positions in the tile set

		// Iterate over the fids used in F_c to make   
		// the candidate connected component.
		submesh.setIdentity(candidate.size(),3);
		int f=0;
		for(auto it2=candidate.begin(); it2!=candidate.end(); it2++)
		{
			for(int i=0; i<3; i++)
			{
				// Populate submesh as a matrix with
				// rows being vids of V that form the tiles
				// of the candidate.
				submesh(f, i) = F_c(*it2,i);

				// Populate V_i as a set of unique vids 
				// encountered while iterating over
				// the tiles which form the candidate.
				if( std::find(V_i.begin(), V_i.end(), F_c(*it2,i)) == V_i.end() )
				{
					V_i.emplace_back(F_c(*it2,i));
				}
			}
			f++;
		}

		// At this point:
		// V_i has the vids of V in the candidate.
		// submesh has #tilesincandidate rows
		// where each row has the 3 vids from 
		// V which make up that tile
		assert(submesh.rows()==6);

		// We need num verts and num edges in candidate
		// connected component for the second test
		std::map<std::pair<int,int>, std::vector<int>> incident_tiles;
		edge_incident_faces(submesh, incident_tiles);

		// Second test
		if(V_i.size()+incident_tiles.size()==V.rows())
		{
			// Create a matrix of vertex positions
			// with new index names
			submesh_vertices.setIdentity(V_i.size(), 3);
			int v=0;
			std::map<int,int> vert_translator;
			for(auto it2=V_i.begin(); it2!=V_i.end(); it2++)
			{
				for(int i=0; i<3; i++)
				{
					vert_translator[*it2] = v;
					submesh_vertices(v,i) = V(*it2,i);
				}
				v++;
			}

			// Rename the vids in submesh to point to 
			// vids in submesh_vertices as opposed 
			// to V, cause they used to have the
			// corners of the current candidate tile.
			for(int f=0; f<submesh.rows(); f++)
			{
				submesh(f,0) = vert_translator[submesh(f,0)];
				submesh(f,1) = vert_translator[submesh(f,1)];
				submesh(f,2) = vert_translator[submesh(f,2)];
			}

			// Subdivide the candidate
			igl::upsample( Eigen::MatrixXd(
				Eigen::MatrixXd(submesh_vertices)), 
				Eigen::MatrixXi(submesh), 
				submesh_vertices, 
				submesh);
			assert(V.rows()==submesh_vertices.rows());
			assert(V.cols()==submesh_vertices.cols());

			// Make sure that every vert in the subdivided
			// candidate can be found in the original V
			int found = true;
			int temp = false;
			std::map<int,bool> ov_visited;
			for(int v=0;v<V.rows();v++)
			{
				ov_visited[v] = false;
			}
			for(int v=0; v<submesh_vertices.rows(); v++)
			{
				for(int ov=0; ov<V.rows();ov++)
				{
					if(
						submesh_vertices(v,0)==V(ov,0) &&
						submesh_vertices(v,1)==V(ov,1) &&
						submesh_vertices(v,2)==V(ov,2) &&
						ov_visited[ov] == false
					){
						temp = true;
						ov_visited[ov] = true;
						std::cout << "Vertex " << v 
											<< " was found at location " 
											<< ov << " !" << std::endl;
					}
				}

				if(temp==true)
				{
					temp = false;
				} else {
					std::cout << "Candidtate failed" << std::endl;
					found = false;
					break;
				}
			}

			// Third (final) test
			if(found) { std::cout << "Gagnant!" << std::endl; }
		} else{ std::cout << "Second test failed." << std::endl; }
	}
};

void sort3(int arr[]) 
{ 
	// Insert arr[1] 
	if (arr[1] < arr[0]) 
		std::swap(arr[0], arr[1]); 

	// Insert arr[2] 
	if (arr[2] < arr[1]) 
	{ 
		std::swap(arr[1], arr[2]); 
		if (arr[1] < arr[0]) 
			std::swap(arr[1], arr[0]); 
	} 
};