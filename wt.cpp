#ifndef IGL_NO_EIGEN
#  include <Eigen/Core>
#endif
#include <iostream>
#include <vector>

using namespace std;

void printstuff()
{
    std::cout << "yoooooo" << std::endl;
};

template <typename DerivedF>
bool edge_incident_faces(
	const Eigen::MatrixBase<DerivedF>& F,
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
	return true;
};

template <typename DerivedF>
bool covering_mesh(
	const Eigen::MatrixBase<DerivedF>& F,
 				Eigen::MatrixBase<DerivedF>& FC
){
	std::map<std::pair<int,int>, std::vector<int>> incident_faces;
  edge_incident_faces(F, incident_faces);
  std::cout << "number of edges: " << incident_faces.size() << std::endl;

	vector< tuple<int, int, int> > tiles;
	for(int f=0; f<F.rows(); f++)
	{
		int v1 = F(f,0);
		int v2 = F(f,1);
		int v3 = F(f,2);

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

			// Neighbouring triangles to current triangle
			int nt1 = incident_faces[e1][0]==f ? incident_faces[e1][1] : f;
			int nt2 = incident_faces[e2][0]==f ? incident_faces[e2][1] : f;
			int nt3 = incident_faces[e3][0]==f ? incident_faces[e3][1] : f;

			// Get tile vertices
			int v1p, v2p, v3p;
			for(int i=0; i<3; i++)
			{
				if(F(nt1, i)!=e1v1 && F(nt1, i)!=e1v2) v1p = F(nt1, i);
				if(F(nt2, i)!=e2v1 && F(nt2, i)!=e2v2) v1p = F(nt2, i);
				if(F(nt3, i)!=e3v1 && F(nt3, i)!=e3v2) v1p = F(nt3, i);
			}
			tiles.push_back(tuple<int, int, int>(v1p,v2p,v3p));
		}

		FC.setIdentity(tiles.size(),3);
		for(int t=0; t < tiles.size(); t++)
		{
			FC(t, 0) = get<0>(tiles[t]);
			FC(t, 1) = get<1>(tiles[t]);
			FC(t, 2) = get<2>(tiles[t]);
		}
	}
	std::cout << "Number of tiles: " << FC.rows() << std::endl;
}