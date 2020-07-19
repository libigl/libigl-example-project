#include "taubin.h"
#include <igl/upsample.h>
#include <iostream>
#include <vector>
#include <map>

// Adapted from Gabriel Taubin's Loop Inverse Subdivsion Algorithm
// Paper: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.21.1125&rep=rep1&type=pdf

// Detect whether the input mesh has subdivision connectivity.
// If so, partitions the input mesh into the coarse 
// mesh and the new vertices after subdividing.
// Returns true if has subdivision connectivity,
// false if not.
bool is_quadrisection(
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXd& V_in,
	Eigen::MatrixXi& F_old,
	Eigen::MatrixXd& V_old,
	Eigen::MatrixXi& F_new,
	Eigen::MatrixXd& V_new
){
  std::cout << "Num verts in input mesh: " << V_in.rows() << std::endl;
  std::cout << "Num faces in input mesh: " << F_in.rows() << std::endl;
  
	// Construct a bunch of tiles
	std::cout << "Begin covering mesh" << std::endl;
  Eigen::MatrixXi tiles; // A tile soup!
  Eigen::MatrixXi covered_faces; // row i is tile i from tiles to 4 fids from F_in
  covering_mesh(F_in, tiles, covered_faces);
	std::cout << "Completed covering mesh" << std::endl;
	// Generate sets of tiles that connect
	std::cout << "Begin connected components" << std::endl;
	std::vector<std::vector<int>> sub_meshes; // Sets of tile soups!
  connected_components(tiles, sub_meshes);
	std::cout << "Completed connected components" << std::endl;
	std::cout << "Number of candidates: " << sub_meshes.size() << std::endl;

	// Find a candidate connected component
	for(auto it=sub_meshes.begin(); it!=sub_meshes.end(); it++)
	{
		// Try to construct a bijection from current
		// connect component to the input mesh
		is_equivalence( F_in, V_in, *it, tiles, covered_faces, F_old, V_old, F_new, V_new );
		// Stop early as soon as we find a
		// successful candidate
		if(F_old.rows()>0 && V_old.rows()>0) 
			return true;
	}
	std::cout << "Mesh does not have subdivision connectivity" << std::endl;
	return false;
};

// Creates tiles from the input mesh,
// And the tile sets covered by each
void covering_mesh(
	const Eigen::MatrixXi& F_in,
 	Eigen::MatrixXi& tiles, // num_tilesx3 matrix of vert ids in V_in
 	Eigen::MatrixXi& covered_faces // # tiles x 4 matrix of fids in F_in
){
	std::map<std::pair<int,int>, std::vector<int>> incident_faces;
  edge_incident_faces(F_in, incident_faces);
  std::cout << "Num edges: " << incident_faces.size() << std::endl;

	std::vector<std::tuple<int, int, int>> found_tiles;
	std::vector<std::tuple<int, int, int, int>> neighbouring_faces;
	for(int f=0; f<F_in.rows(); f++)
	{
		int v1 = F_in(f,0);
		int v2 = F_in(f,1);
		int v3 = F_in(f,2);

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
				if(F_in(nt1, i)!=e1v1 && F_in(nt1, i)!=e1v2) v1p = F_in(nt1, i);
				if(F_in(nt2, i)!=e2v1 && F_in(nt2, i)!=e2v2) v2p = F_in(nt2, i);
				if(F_in(nt3, i)!=e3v1 && F_in(nt3, i)!=e3v2) v3p = F_in(nt3, i);
			}
			found_tiles.push_back(std::tuple<int, int, int>(v1p,v2p,v3p));
			neighbouring_faces.push_back(std::tuple<int, int, int, int>(f,nt1,nt2,nt3));
		}

		tiles.setIdentity(found_tiles.size(),3);
		covered_faces.setIdentity(found_tiles.size(),4);
		for(int t=0; t < found_tiles.size(); t++)
		{
			tiles(t, 0) = std::get<0>(found_tiles[t]);
			tiles(t, 1) = std::get<1>(found_tiles[t]);
			tiles(t, 2) = std::get<2>(found_tiles[t]);

			covered_faces(t, 0) = std::get<0>(neighbouring_faces[t]);
			covered_faces(t, 1) = std::get<1>(neighbouring_faces[t]);
			covered_faces(t, 2) = std::get<2>(neighbouring_faces[t]);
			covered_faces(t, 3) = std::get<3>(neighbouring_faces[t]);
		}
	}
	std::cout << "Num tiles: " << tiles.rows() << std::endl;
};

// Generate all possible connected 
// components from the tiles
void connected_components(
	const Eigen::MatrixXi& tiles,
 	std::vector<std::vector<int>>& sub_meshes // Vector of vectors containing tile ids in a single connected component found
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
		// assert((it->second).size()<=2);
		if((it->second).size()==2)
		{
			int tid1 = it->second[0];
			int tid2 = it->second[1];
			if(where_are_you[tid1]!=where_are_you[tid2])
			{
				assert(where_are_you[tid1]!=NULL && where_are_you[2]!=NULL);

				// where_are_you[tid1]->resize(where_are_you[tid1]->size()+where_are_you[tid2]->size());
				where_are_you[tid1]->insert(
					where_are_you[tid1]->end(),
					where_are_you[tid2]->begin(),
					where_are_you[tid2]->end()
				);

				for(int i=0; i<where_are_you[tid1]->size(); i++)
				{
					if(where_are_you[tid1]->at(i)!=tid2 && where_are_you[tid1]->at(i)!=tid1)
					{
						where_are_you[where_are_you[tid1]->at(i)] = where_are_you[tid1];
					}
				}

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

// Determine whether input connected component
// has a bijection to the original mesh
void is_equivalence(
	const Eigen::MatrixXi& F_in,
	const Eigen::MatrixXd& V_in,
	// Tile indices that make up the candidate mesh.
	const std::vector<int>& candidate,
	// #tiles by 3 matrix of indices into V_in that
	// make up each tile.
	const Eigen::MatrixXi& tiles,
	const Eigen::MatrixXi& covered_faces,
	Eigen::MatrixXi& F_old,
	Eigen::MatrixXd& V_old,
	Eigen::MatrixXi& F_new,
	Eigen::MatrixXd& V_new
){
	// First test
	if(candidate.size()*4==F_in.rows())
	{
		std::cout << "--Begin analyzing the candidate connected component--" << std::endl;

		// Turn candidate into the proper F, V 
		// matrix format so that we can quadrisect it
		Eigen::MatrixXi submesh; // "Faces" from canadidate
		std::vector<int> V_i; // All the vids from V_in that are covered by the candidate
		submesh.setIdentity(candidate.size(),3);
		int t=0;
		for(auto it2=candidate.begin(); it2!=candidate.end(); it2++)
		{
			for(int i=0; i<3; i++)
			{
				// Populate submesh as a matrix with
				// rows being vids of V that form the tiles
				// of the candidate.
				submesh(t, i) = tiles(*it2,i);

				// Populate V_i as a set of unique vids 
				// encountered while iterating over
				// the tiles which form the candidate.
				if( std::find(V_i.begin(), V_i.end(), tiles(*it2,i)) == V_i.end() )
				{
					V_i.emplace_back(tiles(*it2,i));
				}
			}
			t++;
		}

		// At this point:
		// V_i has the vids of V covered by candidate.
		// submesh has #tilesincandidate rows
		// where each row has the 3 vids from 
		// V which make up that tile
		// assert(submesh.rows()==6);

		// We need num verts and num edges in candidate
		// connected component for the second test
		std::map<std::pair<int,int>, std::vector<int>> incident_tiles;
		edge_incident_faces(submesh, incident_tiles);

		// Second test
		if(V_i.size()+incident_tiles.size()==V_in.rows())
		{
			// Make sure that every vert in the subdivided
			// candidate can be found in the original V
			int found = true;
			int temp = false;
			int visited_cout = 0;
			std::map<int,bool> f_visited;
			for(int f=0;f<F_in.rows();f++)
			{
				f_visited[f] = false;
			}

			for(auto it2=candidate.begin(); it2!=candidate.end(); it2++)
			{
				
				for(int i=0; i<4; i++)
				{
					if(!f_visited[covered_faces(*it2,i)])
					{
						f_visited[covered_faces(*it2,i)] = true;
						visited_cout++;
						temp = true;
					}
				}
				if(temp==true)
				{
					temp = false;
				} else {
					std::cout << "Candidate failed" << std::endl;
					found = false;
					break;
				}
			}

			// Third (final) test
			if(found && visited_cout==F_in.rows()) 
			{ 
				std::cout << "Gagnant!" << std::endl; 
				// Create a matrix of vertex positions
				// with new index names
				Eigen::MatrixXd submesh_vertices; // "Vertex positions" from candidate
				submesh_vertices.setIdentity(V_i.size(), 3);
				int v=0;
				std::map<int,int> vert_translator;
				for(auto it2=V_i.begin(); it2!=V_i.end(); it2++)
				{
					for(int i=0; i<3; i++)
					{
						vert_translator[*it2] = v;
						submesh_vertices(v,i) = V_in(*it2,i);
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

				Eigen::MatrixXi F_pre_subdiv = Eigen::MatrixXi(submesh);
				Eigen::MatrixXd V_pre_subdiv = Eigen::MatrixXd(submesh_vertices);

				// Subdivide the candidate
				igl::upsample( Eigen::MatrixXd(
					Eigen::MatrixXd(submesh_vertices)), 
					Eigen::MatrixXi(submesh), 
					submesh_vertices, 
					submesh);
				assert(V_in.rows()==submesh_vertices.rows());
				assert(V_in.cols()==submesh_vertices.cols());

				F_old = F_pre_subdiv;
				V_old = V_pre_subdiv;
				F_new = submesh.block(F_old.rows()-1, 0, submesh.rows()-F_old.rows(), 3);
				V_new = submesh_vertices.block(V_old.rows()-1, 0, submesh_vertices.rows()-V_old.rows(), 3);
				assert(F_new.rows()+F_old.rows() == submesh.rows());
				assert(V_new.rows()+V_old.rows() == submesh_vertices.rows());
			}
		} 
		else
		{ 
			std::cout << "Second test failed." << std::endl; 
		}
	}
};

// Given input mesh, return a map with keys.
// Each key is an ordered pair of vertex 
// ids (v1 comes before v2 if |v1| > |v2|)
// and maps to a list of face ids of 
// all the faces incident to the edge
// formed by v1 and v2.
void edge_incident_faces(
	const Eigen::MatrixXi& F,
	std::map<std::pair<int,int>, std::vector<int>>& incident_faces
){
	for(int f=0; f<F.rows(); f++)
	{
		int v1 = F(f,0);
		int v2 = F(f,1);
		int v3 = F(f,2);

		incident_faces[std::make_pair(std::min(v1,v2),std::max(v1,v2))].push_back(f);
		incident_faces[std::make_pair(std::min(v2,v3),std::max(v2,v3))].push_back(f);
		incident_faces[std::make_pair(std::min(v3,v1),std::max(v3,v1))].push_back(f);
	}
};

// Sort an array of three entries
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