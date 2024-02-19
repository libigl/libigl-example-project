#include <string>
#include <igl/colon.h>
#include <igl/harmonic.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <algorithm>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
// #include <igl/normalize_row_sums.h>
#include <CGAL/Polyhedron_3.h>
#include <igl/vertex_triangle_adjacency.h>
// #include <igl/point_mesh_squared_distance.h>
#include <igl/vertex_triangle_adjacency.h>

#include <igl/is_edge_manifold.h>
#include <CGAL/centroid.h>
#include <igl/loop.h>
#include <igl/writePLY.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <igl/edge_flaps.h>
#include <igl/fast_find_self_intersections.h>
#include <igl/upsample.h>




#include <iostream>
#include <thread>
#include <chrono>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/point_generators_3.h>
#include <igl/collapse_edge.h>
#include <igl/edge_lengths.h>
#include <igl/remove_unreferenced.h>
#include "remesh_botsch.h"

// Define the kernel and triangulation types

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Delaunay_triangulation_3<Kernel> Delaunay;
typedef CGAL::Tetrahedron_3<Kernel> Tetrahedron;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;


// Define a point type
typedef Kernel::Point_3 Point;

double bc_frac = 0.0;
double bc_dir = -0.03;
bool deformation_field = true;
Eigen::MatrixXd V, U, V_bc, U_bc,N,prevU;
Eigen::VectorXd Z;
Eigen::MatrixXi F,prevF,E;
Eigen::VectorXi b;
// Eigen::VectorXi b;
int k=2;
bool grow=true;
bool skip=false;
// Fixed points
std::vector<int> mesh_fixed;
std::vector<int> fixed_points ; // Example: Fixing vertices 0 and 4
std::vector<int> curve_fixed,black_indices;
Eigen::MatrixXd lowPolyVertices;
Eigen::MatrixXd fixedPoints;
double threshold;
int looping=0;

//a function for remeshing the mesh without removing the fixed points and given the threshold
void remesh(double threshold){
	Eigen::VectorXi feature =Eigen::Map<Eigen::VectorXi>(mesh_fixed.data(), mesh_fixed.size());
    Eigen::VectorXd target = Eigen::VectorXd::Constant(V.rows(),threshold);
    remesh_botsch(V,F,target,10,feature,false);
}


// Function to check intersection between a line segment and a triangle
bool lineTriangleIntersection(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
    Eigen::Vector3d e1 = v1 - v0;
    Eigen::Vector3d e2 = v2 - v0;
    Eigen::Vector3d h = p0.cross(p1);
    double a = e1.dot(h);
    
    // Check if the line is parallel to the triangle
    if (a > -1e-6 && a < 1e-6)
        return false;

    double f = 1.0 / a;
    Eigen::Vector3d s = p0 - v0;
    double u = f * s.dot(h);
    
    if (u < 0.0 || u > 1.0)
        return false;

    Eigen::Vector3d q = s.cross(e1);
    double v = f * p1.dot(q);
    
    if (v < 0.0 || u + v > 1.0)
        return false;

    double t = f * e2.dot(q);
    
    if (t > 1e-6)  // Make sure t is positive and not too close to zero
        return true;
    
    return false;
}

// Function to check intersection between a line segment and a mesh
bool lineMeshIntersection( const Eigen::Vector3d& p0, const Eigen::Vector3d& p1) {
    // Iterate through each triangle in the mesh
    for (int i = 0; i < F.rows(); ++i) {
        Eigen::Vector3d x = V.row(F(i, 0));
        Eigen::Vector3d y = V.row(F(i, 1));
        Eigen::Vector3d z = V.row(F(i, 2));

        // Check for intersection with the current triangle
        if (lineTriangleIntersection(p0, p1, x, y, z))
            return true;
    }
    
    // No intersection found with any triangle
    return false;
}

std::vector<int>create_edge_triangle_mapping2( double threshol) {
    std::vector<int> edge_to_triangles;
    // add elements to edge_to_triangles in ascending order of edge_length
    std::vector<std::tuple<double, int>> dis;
    std::cout<<"mesh_fixed"<<mesh_fixed.size()<<std::endl;
    for (int e = 0; e < E.rows(); ++e) {
        
            double edge_length = (V.row(E(e,0)) - V.row(E(e,1))).norm();
            // std::cout<<edge_length<<"--"<<threshold<<std::endl;
            if (edge_length < threshol && std::find(mesh_fixed.begin(), mesh_fixed.end(), E(e,0)) == mesh_fixed.end() &&  std::find(mesh_fixed.begin(), mesh_fixed.end(), E(e,1)) == mesh_fixed.end() ) {
                dis.push_back(std::make_tuple(edge_length, e));
            }
        }
    std::sort(dis.begin(), dis.end());
    
    for(int i=0;i<dis.size();i++){
        edge_to_triangles.push_back(std::get<1>(dis[i]));
    }

    return edge_to_triangles;
    }



void remove_invalid_faces(std::vector<int> invalid_indices)
{   
    Eigen::MatrixXi F_clean;
    

    // Remove faces with invalid indices
    int num_invalid = invalid_indices.size();
    int num_faces = F.rows();
    int num_valid_faces = num_faces - num_invalid;
    F_clean.resize(num_valid_faces, 3);
    int valid_face_idx = 0;

    for (int i = 0; i < num_faces; ++i)
    {
        if (std::find(invalid_indices.begin(), invalid_indices.end(), i) == invalid_indices.end())
        {
            // This face is valid, copy it to the cleaned face matrix
            F_clean.row(valid_face_idx++) = F.row(i);
        }
    }
    F=F_clean;
}


void collapse(double threshol){
using namespace Eigen;
VectorXi EMAP;
  MatrixXi EF,EI;
  
  // Function to reset original mesh and data structures

    igl::edge_flaps(F,E,EMAP,EF,EI);

    std::vector<int> collapse_faces;
    std::vector<int> edges=create_edge_triangle_mapping2(threshol);
     for (int i = 0; i < edges.size(); ++i)
    {
    int e1,e2,f1,f2;
    int e=edges[i];
    if(E(e,0)!=0 && E(e,1)!=0)
    {
    Eigen::RowVector3d  p=(V.row(E(e,0))+V.row(E(e,1)))/2;
    if(igl::collapse_edge(
      e,p,V,F,E,EMAP,EF,EI,e1,e2,f1,f2)){
        collapse_faces.push_back(f1);
        collapse_faces.push_back(f2);
      }
    
    }

    }
    // std::cout<<E.row(e1)<<std::endl<<E.row(e2);

      Eigen::VectorXi I;
    Eigen::MatrixXi newF;
    Eigen::MatrixXd newV;
    igl::remove_unreferenced(V, F, newV, newF, I);
    F=newF;
    V=newV;
    // remove_invalid_faces(collapse_faces);



}




std::vector<int> findInvisiblePoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::RowVector3d& viewpoint)
{
    // Calculate per-vertex normals
    Eigen::MatrixXd VN;
    igl::per_vertex_normals(V, F, VN);

    Eigen::MatrixXd dirVectors = V.rowwise() - viewpoint;

    std::vector<int> invisiblePoints;

    for(int i=0;i<VN.rows();++i){

      if (dirVectors.row(i).dot(VN.row(i))>0.0)
        invisiblePoints.push_back(i);

    }

    return invisiblePoints;
}

int findRowVectorIndex(const Eigen::MatrixXd& matrix, const Eigen::RowVector3d& targetVector)
{
    // Check if matrix has 3 columns (assuming RowVector3d)
    assert(matrix.cols() == 3);

    // Compare each row with the target vector
    for (int i = 0; i < matrix.rows(); ++i)
    {
        if (matrix.row(i) == targetVector)
        {
            return i; // Return the index if found
        }
    }

    return -1; // Return -1 if the vector is not found
}

double compute_largest_edge_length() {


    double max_edge_length = 0.0;

    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1_idx = F(i, j);
            int v2_idx = F(i, (j + 1) % 3);

            const Eigen::Vector3d& v1 = V.row(v1_idx);
            const Eigen::Vector3d& v2 = V.row(v2_idx);

            double edge_length = (v1 - v2).norm();
            max_edge_length = std::max(max_edge_length, edge_length);
        }
    }

    // std::cout<<max_edge_length<<std::endl;
    return max_edge_length;
}

bool isblack(const std::vector<int> black_indices, int index) {
    for (int i = 0; i < black_indices.size(); ++i) {
        if (index == black_indices[i]) {
            return true;
        }
    }
    return false;
}
bool check_visibility_from_vertex(const Eigen::Vector3d& viewer, const Eigen::Vector3d& target, const Eigen::Vector3d& vertex_normal) {
    // Compute the vector from viewer to target
    Eigen::Vector3d view_direction = (target - viewer).normalized();
    // std::cout<<view_direction<<std::endl;
    // Check if the dot product between the view direction and vertex normal is positive
    // If it's negative, the target is not visible from the viewer
    return (view_direction.dot(vertex_normal) < 0.0);
}

std::pair<int, int> find_shortest_edge_between_arrays( Eigen::MatrixXd vertices1,  Eigen::MatrixXd vertices2,  Eigen::MatrixXd& vertex_normal,std::vector<int>curve_fixed,std::vector<int>mesh_fixed) {
    
	int num_vertices1 = vertices1.rows();
    int num_vertices2 = vertices2.rows();

    Eigen::MatrixXd distances(num_vertices1, num_vertices2);
    
    for (int i = 0; i < num_vertices1; ++i) {
        for (int j = 0; j < num_vertices2; ++j) {
            distances(i, j) = (vertices1.row(i) - vertices2.row(j)).norm();
        }
    }
    
    std::vector<std::tuple<double, int, int>> dis;
    for (int i = 0; i < num_vertices1; ++i) {
        for (int j = 0; j < num_vertices2; ++j) {
            if (std::find(mesh_fixed.begin(), mesh_fixed.end(), i) != mesh_fixed.end() ||
                std::find(curve_fixed.begin(), curve_fixed.end(), j) != curve_fixed.end()) {
                // std::cout<<i<<"-"<<j<<std::endl;
                continue; // Skip this iteration
            }
            dis.emplace_back(distances(i, j), i, j);
        }
    }
    std::sort(dis.begin(), dis.end());
    
    for (const auto& d : dis) {
        // std::cout<<vertices1.size()<<"-"<<vertices2.size()<<"-"<<std::get<1>(d)<<"-"<<std::get<2>(d)<<"-"<<vertex_normal.size()<<std::endl;
        // std::cout<<vertices1.row(std::get<1>(d))<<"-"<<vertices2.row(std::get<2>(d))<<"-"<<std::get<1>(d)<<"-"<<std::get<2>(d)<<"-"<<vertex_normal.size()<<std::endl;
        if (check_visibility_from_vertex(vertices2.row(std::get<2>(d)),vertices1.row(std::get<1>(d)), vertex_normal.row(std::get<1>(d)))) {
            Eigen::MatrixXd temp1=vertices1;
            temp1.row(std::get<1>(d))=vertices2.row(std::get<2>(d));
            Eigen::VectorXi EI;
            Eigen::VectorXd p0=vertices1.row(std::get<1>(d));
            Eigen::VectorXd p1=vertices2.row(std::get<2>(d));
            if (!lineMeshIntersection(p0, p1)) {
                return std::make_pair(std::get<1>(d), std::get<2>(d));
            }
            // igl::fast_find_self_intersections(U,F,EI);
            // if(!EI.isZero()){
            //     continue;
            // }
			// if (isblack(black_indices,std::get<2>(d))) {
        
			// 	continue;
			// }
			
            //return std::make_pair(std::get<1>(d), std::get<2>(d));
        }

    }

    std::cout<<" no points"<<std::endl;

    // If no suitable edge is found, return a default value or handle accordingly
    return std::make_pair(-1, -1);
}

void createSphere(const double radius, const Eigen::Vector3d& center, const int num_points=10)
{
    // Initialize vertex and face matrices
    V.resize(num_points * num_points, 3);
    F.resize((num_points - 1) * (num_points - 1) * 2, 3);

    // Generate sphere vertices
    for (int i = 0; i < num_points; ++i)
    {
        for (int j = 0; j < num_points; ++j)
        {
            double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(num_points - 1);
            double phi = M_PI * static_cast<double>(j) / static_cast<double>(num_points - 1);

            V.row(i * num_points + j) << center[0] + radius * sin(phi) * cos(theta),
                                         center[1] + radius * sin(phi) * sin(theta),
                                         center[2] + radius * cos(phi);
        }
    }

    // Generate sphere faces
    int idx = 0;
    for (int i = 0; i < num_points - 1; ++i)
    {
        for (int j = 0; j < num_points - 1; ++j)
        {
            F.row(idx++) << i * num_points + j, (i + 1) * num_points + j, i * num_points + j + 1;
            F.row(idx++) << (i + 1) * num_points + j, (i + 1) * num_points + j + 1, i * num_points + j + 1;
        }
    }
}


double calculateCircumradius(const Point& circumcenter, const Point& vertex) {
    return CGAL::sqrt(CGAL::squared_distance(circumcenter, vertex));
}


bool isPointInsideTetrahedron(const Point& A, const Point& B, const Point& C, const Point& D, const Point& P)
{
    Tetrahedron tetrahedron(A, B, C, D);
    return tetrahedron.has_on_positive_side(P);
}

typedef Eigen::Vector3d Vector3d;

std::pair<double, Vector3d> findSphereThroughPoints(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& D) {
    Vector3d AB = B - A;
    Vector3d AC = C - A;
    Vector3d AD = D - A;

    // Calculate the 3x3 matrix determinant
    double detA = AB.cross(AC).dot(AD);

    // Calculate the coefficients of the quadratic equation
    double a = 2.0 * (AB.dot(AC.cross(AD)));
    double b = AB.squaredNorm() * (AD.cross(AC)).dot(AB) +
               AC.squaredNorm() * (AB.cross(AD)).dot(AC) +
               AD.squaredNorm() * (AC.cross(AB)).dot(AD);
    double c = AB.cross(AC).squaredNorm() * AD.squaredNorm() +
               AC.cross(AD).squaredNorm() * AB.squaredNorm() +
               AD.cross(AB).squaredNorm() * AC.squaredNorm();

    // Calculate the radius of the sphere
    double radius = std::sqrt(b * b - 4.0 * a * c) / (2.0 * std::fabs(detA));

    // Calculate the center of the sphere
    Vector3d center = A + (AB.cross(AC).cross(AD) * (b / (2.0 * std::fabs(detA))));

    return std::make_pair(radius, center);
}
std::pair<double, Eigen::Vector3d> computesphere(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& D) {
    Vector3d BA = A - B;
    Vector3d BC = C - B;
    Vector3d BD = D - B;

    double dot_BA_BC = BA.dot(BC);
    double dot_BD_BC = BD.dot(BC);

    Vector3d BC_cross_BA = BC.cross(BA);
    Vector3d BD_cross_BC = BD.cross(BC);

    // Calculate the circumradius
    double circumradius = 0.5 * BC.norm() * BC_cross_BA.norm() / BC_cross_BA.cross(BD_cross_BC).norm();

    // Calculate the circumcenter
    Vector3d circumcenter = B + dot_BD_BC * BC_cross_BA / BC_cross_BA.squaredNorm();

    return std::make_pair(circumradius, circumcenter);
}
void transformIcosahedron(const std::string& objFilename, double desiredRadius, const Eigen::Vector3d& center) {
    Eigen::MatrixXd originalV;
    Eigen::MatrixXi originalF;

    // Read the icosahedron from the OBJ file
    igl::readOBJ(objFilename, originalV, originalF);
    Eigen::Vector3d origin(0.0,0.0,0.0);
    // Calculate the original size of the icosahedron
    
    double originalSize = originalV.row(0).norm();

    // Scale the vertices to achieve the desired radius based on the original size
    originalV *= (desiredRadius / originalSize);

    // Translate the vertices to the desired center
    originalV.rowwise() += center.transpose();

    V = originalV;
    F = originalF;
}

bool isOriginInsideSphere(const Point& sphereCenter, double sphereRadius) {
    Point origin(0.0, 0.0, 0.0);

    double distanceSquared = CGAL::squared_distance(origin, sphereCenter);
    double radiusSquared = sphereRadius * sphereRadius;

    return distanceSquared < radiusSquared;
}

void computeCircumradiusAndCenter(const Eigen::MatrixXd& points1) {
    // Eigen::MatrixXi DT;
    
    Point largecent ;
    std::vector<Point> points;
    double largestCircumradius = -1.0;
    Eigen::Vector3d largestCircumcenter(0.0, 0.0, 0.0);
    
    for (int i = 0; i < points1.rows(); ++i) {
        Point p(points1(i, 0), points1(i, 1), points1(i, 2));
        points.push_back(p);
    }

    Delaunay dt(points.begin(), points.end());

    // Iterate through the tetrahedra of the Delaunay triangulation
    for (auto it = dt.finite_cells_begin(); it != dt.finite_cells_end(); ++it) {
        // Access the vertices of each tetrahedron
        const Point& d1 = it->vertex(0)->point();
        const Point& d2 = it->vertex(1)->point();
        const Point& d3 = it->vertex(2)->point();
        const Point& d4 = it->vertex(3)->point();

        // std::cout<<it->vertex(0)->point()<<std::endl;

        Eigen::MatrixXd tetrahedronVertices(4, 3);

        tetrahedronVertices.row(0) = Eigen::Vector3d(d1.x(),d1.y(),d1.z());
        tetrahedronVertices.row(1) = Eigen::Vector3d(d2.x(),d2.y(),d2.z());
        tetrahedronVertices.row(2) = Eigen::Vector3d(d3.x(),d3.y(),d3.z());
        tetrahedronVertices.row(3) = Eigen::Vector3d(d4.x(),d4.y(),d4.z());

        // std::cout<<tetrahedronVertices.row(0)<<std::endl;

        // Eigen::Vector3d origin(0.0,0.0,0.0);
        Point p(0.0,0.0,0.0);

        if(isPointInsideTetrahedron(d1, d2, d3, d4,p)) {
                // std::pair<double, Eigen::Vector3d> sphere=computesphere(tetrahedronVertices.row(0),tetrahedronVertices.row(1),tetrahedronVertices.row(2),tetrahedronVertices.row(3));
                // addTetrahedronToMesh(tetrahedronVertices);
                // int numVertices = V.rows();

                // // Add tetrahedron vertices to V
                // V.conservativeResize(numVertices + 4, 3);
                // V.row(numVertices) = Eigen::Vector3d(d1.x(), d1.y(), d1.z());
                // V.row(numVertices + 1) = Eigen::Vector3d(d2.x(), d2.y(), d2.z());
                // V.row(numVertices + 2) = Eigen::Vector3d(d3.x(), d3.y(), d3.z());
                // V.row(numVertices + 3) = Eigen::Vector3d(d4.x(), d4.y(), d4.z());

                // // Add tetrahedron faces to F
                // int numFaces = F.rows();
                // F.conservativeResize(numFaces + 4, 3);
                // F.row(numFaces) << numVertices, numVertices + 1, numVertices + 2;
                // F.row(numFaces + 1) << numVertices, numVertices + 2, numVertices + 3;
                // F.row(numFaces + 2) << numVertices, numVertices + 3, numVertices + 1;
                // F.row(numFaces + 3) << numVertices + 1, numVertices + 2, numVertices + 3;
                // largecent=circumcenter;
                Point circumcenter = CGAL::centroid(d1, d2, d3, d4);
        
                // double circumradius = CGAL::sqrt(CGAL::squared_distance(circumcenter, d1));

                largestCircumcenter=Vector3d(circumcenter.x(),circumcenter.y(),circumcenter.z());
                largestCircumradius=CGAL::sqrt(CGAL::squared_distance(circumcenter, d1));
                // std::cout<<largestCircumcenter;

                // computespahere(tetrahedronVertices.row(0),tetrahedronVertices.row(1),tetrahedronVertices.row(2),tetrahedronVertices.row(3));
        }

        }
    
    // createIcosahedron(largestCircumradius,largestCircumcenter);
    transformIcosahedron("/home/anandhu/Desktop/libigl-example-project/sphere.obj",largestCircumradius, largestCircumcenter);

    // createSphere(largestCircumradius,largestCircumcenter);

}

bool endsWith(const std::string& str, const std::string& suffix) {
    if (str.length() >= suffix.length()) {
        return (str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0);
    } else {
        return false;
    }
}


std::map<std::pair<int, int>, std::vector<int>> create_edge_triangle_mapping( double threshold) {
    std::map<std::pair<int, int>, std::vector<int>> edge_to_triangles;

    for (int triangle_idx = 0; triangle_idx < F.rows(); ++triangle_idx) {
        const Eigen::RowVector3i& triangle = F.row(triangle_idx);

        for (int j = 0; j < 3; ++j) {
            int v0 = triangle(j);
            int v1 = triangle((j + 1) % 3);

            double edge_length = (V.row(v0) - V.row(v1)).norm();
            // std::cout<<edge_length<<"--"<<threshold<<std::endl;
            if (edge_length > threshold) {
                std::pair<int, int> edge = std::make_pair(std::min(v0, v1), std::max(v0, v1));
                edge_to_triangles[edge].push_back(triangle_idx);
            }
        }
    }

    return edge_to_triangles;
}



void subdivide_mesh_with_mapping(double threshol) {
    
    // std::cout<<threshol;
    int count=0;
    std::map<std::pair<int, int>, std::vector<int>> edge_triangle_mapping =
        create_edge_triangle_mapping(threshol);

    while (edge_triangle_mapping.size()>0 )
    {
        Eigen::MatrixXi old_F = F;
    std::cout<<"subdivision STARTED";
    igl::per_vertex_normals(V, F, N);
    std::vector<int> divide_triangle;
        /* code */
    
    
    int new_vertex_index = V.rows();
    F.resize(0, 3);

    for (const auto& entry : edge_triangle_mapping) {
        
        bool subdivide = true;

        const std::pair<int, int>& edge = entry.first;
        const std::vector<int>& triangle_indices = entry.second;

        if (std::find(divide_triangle.begin(), divide_triangle.end(), triangle_indices[0]) != divide_triangle.end()) {
            subdivide = false;
        }

        if (std::find(divide_triangle.begin(), divide_triangle.end(), triangle_indices[1]) != divide_triangle.end()) {
            subdivide = false;
        }

        if (subdivide) {

            const std::pair<int, int>& edge = entry.first;
            const std::vector<int>& triangle_indices = entry.second;
            Eigen::RowVector3d midpoint = (V.row(edge.first) + V.row(edge.second)) / 2.0;
            // std::cout<<midpoint<<std::endl;
            V.conservativeResize(new_vertex_index + 1, Eigen::NoChange);
            V.row(new_vertex_index) = midpoint;
            
            Eigen::Vector3d normal_sum = N.row(edge.first) + N.row(edge.second);
            Eigen::RowVector3d midpoint_normal = normal_sum.normalized();
            N.conservativeResize(new_vertex_index + 1, Eigen::NoChange);
            N.row(new_vertex_index) = midpoint_normal;

        //     std::cout<<"subdivision done";
        // exit(0);

            for (int triangle_idx : triangle_indices) {
                const Eigen::RowVector3i& triangle = old_F.row(triangle_idx);
                divide_triangle.push_back(triangle_idx);
                std::vector<int> non_edge_vertex_idx;
                int index_edge_first,index_edge_second;
                for (int j = 0; j < 3; ++j) {
                    if (triangle[j] != edge.first && triangle[j] != edge.second) {
                        non_edge_vertex_idx.push_back(triangle[j]);
                    }
                    if (triangle[j] == edge.first) {
                        index_edge_first = j;
                    }
                    if (triangle[j] == edge.second) {
                        index_edge_second = j;
                    }
                }
                Eigen::RowVector3i new_tri1,new_tri2;
                new_tri1= triangle;
                new_tri2 =triangle;
                new_tri1[index_edge_second]=new_vertex_index;
                new_tri2[index_edge_first]=new_vertex_index;

                // std::cout<<triangle<<"--"<<new_tri1<<"--"<<new_tri2<<std::endl;
                
                F.conservativeResize(F.rows() + 2, Eigen::NoChange);
                F.row(F.rows() - 2) =new_tri1;
                F.row(F.rows() - 1) =new_tri2;

            }
            new_vertex_index++;
            
        }
    }

    for (int i = 0; i < old_F.rows(); i++) {

        bool triangle_added = false;

        for (int j = 0; j < divide_triangle.size(); j++) {
            if (i == divide_triangle[j]) {
                triangle_added = true;
                break;
            }
        }
        

        if (!triangle_added) {
            F.conservativeResize(F.rows() + 1, Eigen::NoChange);
            F.row(F.rows() - 1) = old_F.row(i);
        }
        
    }

    edge_triangle_mapping =
        create_edge_triangle_mapping(threshol);
    count++;
    }

    
    // U=V;
    
}


// Initialize function
bool initialize(const  std::string file_path,double scale=1) {
    // Load mesh (replace "your_ply" with your mesh file)
    // Load mesh (replace "your_ply" with your mesh file)
    // if (!igl::readOBJ("/home/anandhu/Desktop/libigl-example-project/sphere.obj", V, F)) {
    //     std::cerr << "Error: Couldn't load " << std::endl;
    //     return false;
    // }

    // Eigen::MatrixXd highPolyVertices;
    // Eigen::MatrixXi highPolyFaces;
    // igl::readPLY("/home/anandhu/Desktop/libigl-example-project/chick.ply", highPolyVertices, highPolyFaces);

    // Load the low poly mesh (replace "low_poly_obj" with your low poly mesh file)
    
    Eigen::MatrixXi lowPolyFaces;

    // Check if the file path ends with ".obj"
    if (endsWith(file_path,".obj")) {

    if (!igl::readOBJ(file_path, lowPolyVertices, lowPolyFaces)) {
        std::cout << "Error: Failed to read the OBJ file." << std::endl;
        return 1; // Return an error code indicating file read failure
    }
    }
     else if(endsWith(file_path,".ply")){
        igl::readPLY(file_path, lowPolyVertices, lowPolyFaces);
    }
    else if(endsWith(file_path,".xyz")){

        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file!" << std::endl;
            return -1;
        }

    // Initialize Eigen matrix to store points


    std::string line;
    int count=0;
    while (std::getline(file, line) and count<150) {
        std::istringstream iss(line);
        double x, y, z;
        if (!(iss >> x >> y >> z)) {
            std::cerr << "Error: Invalid file format!" << std::endl;
            return -1;
        }
        lowPolyVertices.conservativeResize(lowPolyVertices.rows() + 1, 3);
        lowPolyVertices.row(lowPolyVertices.rows() - 1) << x, y, z;
        count++;
    }
    file.close();
    }

    

    computeCircumradiusAndCenter(lowPolyVertices);
    // igl::loop(V,F,V,F);
    threshold=compute_largest_edge_length();
    threshold=threshold/scale;
	remesh(threshold);
    // subdivide_mesh_with_mapping(threshold);
    // igl::upsample(V,F,V,F);
	// igl::upsample(V,F,V,F);
	// igl::loop(V,F,V,F);
    // std::cout<<threshold<<"threshold"<<std::endl;
        
    
    return true;
}


std::vector<int> findIndices(const std::vector<int> fixed_indices ) {
    std::vector<int> indices;
    for (int i = 0; i < fixed_indices.size(); ++i) {
        for (int j = 0; j < V.rows(); ++j) {
            if ((V.row(j) - lowPolyVertices.row(fixed_indices[i])).norm() < 1e-6) {
                indices.push_back(j);
                break;  // Move to the next lowPolyVertex
            }
        }
    }
    return indices;
}

std::vector<int> findNeighboringVertices(int vertexIndex) {
    std::unordered_set<int> neighborVertices;

    // Iterate over all the faces
    for (int i = 0; i < F.rows(); ++i) {
        // Check if the vertex is present in the current face
        for (int j = 0; j < F.cols(); ++j) {
            if (F(i, j) == vertexIndex) {
                // Add the other vertices of the face to the set
                for (int k = 0; k < F.cols(); ++k) {
                    if (k != j) {
                        neighborVertices.insert(F(i, k));
                    }
                }
                break; // No need to check other vertices of this face
            }
        }
    }

    // Convert unordered_set to vector
    std::vector<int> neighbors(neighborVertices.begin(), neighborVertices.end());
    return neighbors;
}

bool hasNaN(const Eigen::MatrixXd& G)
{
    // Iterate through all vertices
    for (int i = 0; i < G.rows(); ++i)
    {
        // Check if any coordinate of the vertex is NaN
        if (G.row(i).hasNaN())
        {
            return true; // Found a NaN value
        }
    }
    return false; // No NaN value found
}

bool pre_draw(igl::opengl::glfw::Viewer& viewer)
{   
    if(!grow){ 
        
        viewer.core().is_animating = false;
        return false;}
    
    using namespace Eigen;

    Eigen::MatrixXd normals;
    igl::per_vertex_normals(V, F, normals);
    mesh_fixed=findIndices(curve_fixed);
    std::pair<int, int> result=find_shortest_edge_between_arrays(V,lowPolyVertices,normals,curve_fixed,mesh_fixed);
	std::cout<<(curve_fixed.size())<<std::endl;
  std::cout<<(result.first)<<std::endl;
	std::cout<<(result.second)<<std::endl;
    if(result.first==-1) return false;
    
//   V.row(result.first)= lowPolyVertices.row(result.second);
//   std::cout<<lowPolyVertices.row(result.second);
    // subdivide_mesh_with_mapping(4*threshold);


	// curve_fixed.push_back(result.second);

	mesh_fixed.push_back(result.first);
  std::vector<int> invisible;
  if (curve_fixed.size()<50){
    invisible=findInvisiblePoints(V,F,lowPolyVertices.row(result.second));
  }
//   for (int c = 0; c < mesh_fixed.size()-1; ++c) {
//         int index=findRowVectorIndex(V,lowPolyVertices.row(curve_fixed[c]));
//         mesh_fixed.push_back(index);
//     }

    VectorXi S;
    S.resize(U.rows());

    S.setConstant(-1);
    
    for (int i = 0; i < invisible.size(); i++) {
            S(invisible[i]) = 0;  // Casting double to int, replace with your logic
            // std::cout<<mesh_fixed[i]<<std::endl;
        }
    for (int i = 0; i < mesh_fixed.size()-1; i++) {
            S(mesh_fixed[i]) = 0;  // Casting double to int, replace with your logic
            // std::cout<<mesh_fixed[i]<<std::endl;
    }
    S(result.first) = 1;
    int count_minus_one = 0;
    for (int i = 0; i < S.size(); ++i) {
        if (S[i] == -1) {
            count_minus_one++;
        }
    }
    int count_zero = 0;
    for (int i = 0; i < S.size(); ++i) {
        if (S[i] == 0) {
            count_zero++;
        }
    }

    std::vector<int> neighs= findNeighboringVertices(result.first);
    // std::cout << "Count of elements with fixed neighbour v: " << fixed_neighs << std::endl;
    int fixed_neighs=0;
    for (int i = 0; i < neighs.size(); ++i) {
        if (S[neighs[i]] == -1) {
            fixed_neighs++;
        }
    }

    // Print the count
    // std::cout << "Count of elements with value -1 in vector S: " << count_minus_one << std::endl;
    std::cout << "Count of elements with fixed neighbour v: " << fixed_neighs << std::endl;
    igl::colon<int>(0,V.rows()-1,b);
    
    try{
    b.conservativeResize(std::stable_partition( b.data(), b.data()+b.size(),
    [&S](int i)->bool{return S(i)>=0;})-b.data());
    }
    catch(const std::exception& e) {

        std::cout<<"conservative issue"<<std::endl;
     }
    //  std::cout<<"bsize"<<std::to_string(b.size())<<std::endl;


    // if (curve_fixed.size()>=21 )
    //    {
    //     Eigen::MatrixXd plot(2,3);
    //     plot.row(0)=lowPolyVertices.row(result.second);
    //     plot.row(1)=V.row(result.first);
    //     if(grow)
    //     for(int bi = 0;bi<b.size();bi++){
    //         plot.conservativeResize(plot.rows()+1,3);
    //         plot.row(plot.rows()-1)=V.row(b(bi));
    //     }
    
    //     viewer.data().clear();
    //     viewer.data().set_mesh(U, F);
    //     viewer.data().set_points(lowPolyVertices, Eigen::RowVector3d(1.0, 0.0, 0.0));
    //     viewer.data().set_points(plot, Eigen::RowVector3d(0.0, 1.0, 0.0));
        
    // return false;
    // }
    
    curve_fixed.push_back(result.second);

    // Boundary conditions directly on deformed positions
    U_bc.resize(b.size(),V.cols());
    V_bc.resize(b.size(),V.cols());
    for(int bi = 0;bi<b.size();bi++)
    {
        V_bc.row(bi) = V.row(b(bi));
        switch(S(b(bi)))
        {
        case 0:
            // Don't move handle 0
            U_bc.row(bi) = V.row(b(bi));
            break;
        case 1:
            // Don't move handle 0
            U_bc.row(bi) = lowPolyVertices.row(result.second);;
            break;
        

        }
    }

    // if(viewer.core().is_animating)
    // {
    //     bc_frac += bc_dir;
    //     bc_dir *= (bc_frac>=1.0 || bc_frac<=0.0?-1.0:1.0);
    // }

    
    // Try to perform harmonic deformation
    
    try {
        // Attempt harmonic deformation
        if (igl::harmonic(V, F, b, U_bc, k, U)) {
            // Harmonic deformation succeeded
            std::cout << "Harmonic deformation done" << std::endl;
            Eigen::VectorXi EI;
            Eigen::MatrixXd EV;
            Eigen::MatrixXi IF,EE;
			// igl::fast_find_self_intersections(U,F,EI);
            // if(!EI.isZero())
            // {   std::cout << "harmonic-interssecttt" << std::endl;
            // igl::writeOBJ("/home/anandhu/Desktop/libigl-example-project/output/intersection.obj", U, F);
            // std::ofstream outFile("EI.txt");

			// if(black_indices.size()>10){
			// 	black_indices.clear();
			// 	curve_fixed.pop_back();
			// 	threshold=threshold*0.7;
			// 	// subdivide_mesh_with_mapping(threshold);
            //     // igl::upsample(V,F,V,F);
			// 	U=V;
			// 	return false;
			// }
			// else{
			// 	black_indices.push_back(result.second);
			// 	curve_fixed.pop_back();
			// 	// remesh(threshold);
			// 	U=V;
			// 	return false;

			// }
            // }
			// else{
			// 	black_indices.clear();
			// }
			if(hasNaN(U)){
				std::cout << "harmonic-NAN" << std::endl;
				exit(0);
			}
			
				   
				// igl::writeOBJ("/home/anandhu/Desktop/libigl-example-project/output/outnanU.ply", U_bc);

			
		} else {
			// Harmonic deformation failed, handle the errorthreft
			//
			U = V;
			U.row(result.first) = lowPolyVertices.row(result.second);
			std::cerr << "Harmonic deformation failed" << std::endl;
		}
	} catch (const std::exception& e) {
		// Handle the exception
		U = V;
		U.row(result.first) = lowPolyVertices.row(result.second);
		std::cerr << "Exception during harmonic deformation: " << e.what() << std::endl;
	}
	// if (curve_fixed.size()==22 )
	//    {
	//     Eigen::MatrixXd plot(2,3);
	//     plot.row(0)=lowPolyVertices.row(result.second);
	//     plot.row(1)=V.row(result.first);
	//     if(grow)
	//     for(int bi = 0;bi<b.size();bi++){
	//         plot.conservativeResize(plot.rows()+1,3);
	//         plot.row(plot.rows()-1)=V.row(b(bi));
	//     }
	
	//     std::cout<<lowPolyVertices.row(result.second)<<std::endl;
	//     std::cout<<((V.row(result.first)-lowPolyVertices.row(result.second)).norm()>= 1e-6)<<std::endl;;
	//     viewer.data().clear();
	//     viewer.data().set_mesh(U, F);
	//     viewer.data().set_points(lowPolyVertices, Eigen::RowVector3d(1.0, 0.0, 0.0));
	//     viewer.data().set_points(plot, Eigen::RowVector3d(0.0, 1.0, 0.0));
	//     igl::writeOBJ("/home/anandhu


                        
    // if (curve_fixed.size()==22 )
    //    {
    //     Eigen::MatrixXd plot(2,3);
    //     plot.row(0)=lowPolyVertices.row(result.second);
    //     plot.row(1)=V.row(result.first);
    //     if(grow)
    //     for(int bi = 0;bi<b.size();bi++){
    //         plot.conservativeResize(plot.rows()+1,3);
    //         plot.row(plot.rows()-1)=V.row(b(bi));
    //     }
    
    //     std::cout<<lowPolyVertices.row(result.second)<<std::endl;
    //     std::cout<<((V.row(result.first)-lowPolyVertices.row(result.second)).norm()>= 1e-6)<<std::endl;;
    //     viewer.data().clear();
    //     viewer.data().set_mesh(U, F);
    //     viewer.data().set_points(lowPolyVertices, Eigen::RowVector3d(1.0, 0.0, 0.0));
    //     viewer.data().set_points(plot, Eigen::RowVector3d(0.0, 1.0, 0.0));
    //     igl::writeOBJ("/home/anandhu/Desktop/libigl-example-project/output/outnanV.obj", V, F);
    //     igl::writeOBJ("/home/anandhu/Desktop/libigl-example-project/output/outnanU.obj", U, F);
    //     exit(0);
        
    // return false;
    // }
    



//     fixedPoints.conservativeResize(0, 3);
//     for (int c = 0; c < invisible.size(); ++c) {
//         // mesh_fixed.push_back(invisible[c]);
//         fixedPoints.conservativeResize(fixedPoints.rows() + 1, 3);
//         fixedPoints.row(fixedPoints.rows() - 1)=U.row(invisible[c]);
//     }
//     for (int c = 0; c < mesh_fixed.size(); ++c) {
//         // mesh_fixed.push_back(invisible[c]);
//         fixedPoints.conservativeResize(fixedPoints.rows() + 1, 3);
//         fixedPoints.row(fixedPoints.rows() - 1)=U.row(mesh_fixed[c]);
//     }
    
//     MatrixXd C(F.rows(),3);
//     RowVector3d purple(80.0/255.0,64.0/255.0,255.0/255.0);
//     RowVector3d gold(255.0/255.0,228.0/255.0,58.0/255.0);
//     for(int f = 0;f<F.rows();f++)
//   {
//     if( S(F(f,0))>=0&& S(F(f,1))>=0 && S(F(f,2))>=0)
//     {
//       C.row(f) = purple;
//     }else
//     {
//       C.row(f) = gold;
//     }
//   }
    viewer.data().clear();
    viewer.data().set_mesh(U, F);
    viewer.data().set_points(lowPolyVertices, Eigen::RowVector3d(0.0, 1.0, 0.0));
    // viewer.data().set_colors(C);
    V=U;
    
    // viewer.data().set_points(fixedPoints, Eigen::RowVector3d(1.0, 0.0, 0.0));
    
    // collapse(threshold/2.1);
    // mesh_fixed=findIndices(curve_fixed);
    // collapse(threshold/2.1);
	remesh(threshold);
    Eigen::VectorXi EI;
    Eigen::MatrixXd EV;
    Eigen::MatrixXi IF,EE;
    // igl::fast_find_self_intersections(V,F,EI);
    //         if(!EI.isZero())
    //         {   std::cout << "collapse-interssecttt" << std::endl;
    //         igl::writeOBJ("/home/anandhu/Desktop/libigl-example-project/output/intersection.obj", V, F);
    //         std::ofstream outFile("EI.txt");

    //             exit(0);
    //         }
    // igl::writeOBJ("/home/anandhu/Desktop/libigl-example-project/output/before-out"+std::to_string(curve_fixed.size())+".obj", V, F);

    // subdivide_mesh_with_mapping(threshold);

//    igl::upsample(V,F,V,F);


    if(!igl::is_edge_manifold(F))
    {
        std::cout<<"non-manifoldddd"<<std::endl;
        F.resize(prevF.rows(),3);
        F=prevF;
        V.resize(prevU.size(),3);
        V=prevU;
        exit(0);
    }

    U=V;
    
    
    // viewer.data().compute_normals();
    igl::writeOBJ("/home/anandhu/Desktop/libigl-example-project/output/after-out"+std::to_string(curve_fixed.size())+".obj", V, F);
    // std::chrono::seconds duration(5);
    // std::this_thread::sleep_for(duration);
    return false;
}


// bool pre_draw2(igl::opengl::glfw::Viewer& viewer)
// {
//     if (curve_fixed.size()>=lowPolyVertices.size())
//        {
//     return false;
//     }
//     prevU=V;
//     prevF=F;
//     using namespace Eigen;

//     Eigen::MatrixXd normals;
//     igl::per_vertex_normals(V, F, normals);
//     mesh_fixed=findIndices(curve_fixed);
//     std::pair<int, int> result=find_shortest_edge_between_arrays(V,lowPolyVertices,normals,curve_fixed,mesh_fixed);
// 	std::cout<<(curve_fixed.size())<<std::endl;
//   std::cout<<(result.first)<<std::endl;
// 	std::cout<<(result.second)<<std::endl;
//     if(result.first==-1) return false;


// 	curve_fixed.push_back(result.second);

// 	mesh_fixed.push_back(result.first);
  
// //   U.row(result.first)= lowPolyVertices.row(result.second);

//   std::vector<int> invisible=findInvisiblePoints(V,F,lowPolyVertices.row(result.second));
// //   for (int c = 0; c < mesh_fixed.size()-1; ++c) {
// //         int index=findRowVectorIndex(V,lowPolyVertices.row(curve_fixed[c]));
// //         mesh_fixed.push_back(index);
// //     }

//     VectorXi S;
//     S.resize(U.rows());

//     //   igl::readDMAT(TUTORIAL_SHARED_PATH "/decimated-max-selection.dmat",S);
//     S.setConstant(-1);
//     S(result.first) = 1;
//     for (int i = 0; i < invisible.size(); i++) {
//             S(invisible[i]) = 0;  // Casting double to int, replace with your logic
//             // std::cout<<mesh_fixed[i]<<std::endl;
//         }
//     for (int i = 0; i < mesh_fixed.size()-1; i++) {
//             S(mesh_fixed[i]) = 0;  // Casting double to int, replace with your logic
//             // std::cout<<mesh_fixed[i]<<std::endl;
//         }

//     int count_minus_one = 0;
//     for (int i = 0; i < S.size(); ++i) {
//         if (S[i] == -1) {
//             count_minus_one++;
//         }
//     }

//     igl::colon<int>(0,V.rows()-1,b);
    
//     b.conservativeResize(std::stable_partition( b.data(), b.data()+b.size(),
//     [&S](int i)->bool{return S(i)>=0;})-b.data());

//     // Boundary conditions directly on deformed positions
//     U_bc.resize(b.size(),V.cols());
//     V_bc.resize(b.size(),V.cols());
//     for(int bi = 0;bi<b.size();bi++)
//     {
//         V_bc.row(bi) = V.row(b(bi));
//         switch(S(b(bi)))
//         {
//         case 0:
//             // Don't move handle 0
//             U_bc.row(bi) = V.row(b(bi));
//             break;
//         case 1:
//             // Don't move handle 0
//             U_bc.row(bi) = lowPolyVertices.row(result.second);;
//             break;
        

//         }
//     }
    
//      if(count_minus_one!=0){
//      igl::harmonic(V, F, b, U_bc, k, U);
//     }
//     else{
//         U=V;
//         U.row(result.first)=lowPolyVertices.row(result.second);
//     }
//      fixedPoints.conservativeResize(0, 3);
//     for (int c = 0; c < invisible.size(); ++c) {
//         // mesh_fixed.push_back(invisible[c]);
//         fixedPoints.conservativeResize(fixedPoints.rows() + 1, 3);
//         fixedPoints.row(fixedPoints.rows() - 1)=U.row(invisible[c]);
//     }
//     for (int c = 0; c < mesh_fixed.size(); ++c) {
//         // mesh_fixed.push_back(invisible[c]);
//         fixedPoints.conservativeResize(fixedPoints.rows() + 1, 3);
//         fixedPoints.row(fixedPoints.rows() - 1)=U.row(mesh_fixed[c]);
//     }
    
//     MatrixXd C(F.rows(),3);
//     RowVector3d purple(80.0/255.0,64.0/255.0,255.0/255.0);
//     RowVector3d gold(255.0/255.0,228.0/255.0,58.0/255.0);
//     for(int f = 0;f<F.rows();f++)
//   {
//     if( S(F(f,0))>=0&& S(F(f,1))>=0 && S(F(f,2))>=0)
//     {
//       C.row(f) = purple;
//     }else
//     {
//       C.row(f) = gold;
//     }
//   }
//     viewer.data().clear();
//     viewer.data().set_mesh(U, F);
//     viewer.data().set_points(lowPolyVertices, Eigen::RowVector3d(0.0, 1.0, 0.0));
//     viewer.data().set_colors(C);
//     V=U;
    
//     // viewer.data().set_points(fixedPoints, Eigen::RowVector3d(1.0, 0.0, 0.0));
//     subdivide_mesh_with_mapping(threshold);
//     // edge_collapse_mesh(threshold/2.1);
//     U=V;
//     // viewer.data().compute_normals();
//     igl::writeOBJ("/home/anandhu/Desktop/libigl-example-project/output/out"+std::to_string(curve_fixed.size())+".obj", V, F);
//     // std::chrono::seconds duration(5);
//     // std::this_thread::sleep_for(duration);
//     return false;
// }

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mods)
{
    switch (key)
    {
    case ' ':
        viewer.core().is_animating = !viewer.core().is_animating;
        return true;
    case 'D':
    case 'd':
        deformation_field = !deformation_field;
        return true;
    case '.':
      k++;
      k = (k>4?4:k);
      break;
    case ',':
      k--;
      k = (k<1?1:k);
      break;
    case 'g':
    case 'G':
        grow=true;
//    while (curve_fixed.size()<=lowPolyVertices.size() )
//        {
        // pre_draw2(viewer);
    // }
        break;
    case 's':
    case 'S':
        grow=false;
    // while(grow){
        // inflate(viewer);
    // }
        break;
    case 'b':
    case 'B':
        U=V=prevU;
        F=prevF;
        break;
    case 'p':
    case 'P':
      std::cout<<V.row(20);
      break;
    }
    
    return false;
}

int main(int argc, char* argv[])
{
    using namespace Eigen;
    using namespace std;
    // igl::readOBJ("sphere.obj", V, F);  // Replace "your_mesh.obj" with the path to your mesh file
    // Check if the user provided the file path argument
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <file_path>" << endl;
        return 1; // Return an error code indicating incorrect usage
    }

    char* endptr;
    // Read the file path from the command-line arguments
    const std::string file_path = argv[1];
    double scale = std::strtod(argv[2], &endptr);

    initialize(file_path,scale);
    igl::per_vertex_normals(V, F, N);
    U = V;

    // Plot the mesh with pseudocolors
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(U, F);
    viewer.data().set_points(lowPolyVertices, Eigen::RowVector3d(1.0, 0.0, 0.0)); // Points color: Red

    viewer.data().show_lines = false;
    viewer.core().viewport = Eigen::Vector4f(0, 0, 1800, 1600);
    // viewer.data().set_colors(C);
    viewer.core().trackball_angle = Eigen::Quaternionf(sqrt(2.0), 0, sqrt(2.0), 0);
    viewer.core().trackball_angle.normalize();
    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &key_down;
    viewer.core().is_animating = true;
    viewer.core().animation_max_fps = 30.;
    
    cout <<
        "Press [space] to toggle deformation." << endl <<
        "Press 'd' to toggle between biharmonic surface or displacements." << endl;
    viewer.launch();
    //  while (!viewer.core().is_animating) {
    //     // Draw the current frame
    //     viewer.draw();
    //     viewer.core().is_animating=false;
    // }

}