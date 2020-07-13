#include <igl/PI.h>
#include <igl/upsample.h>
#include <iostream>
#include <vector>
#include <map>

// Wei Shengyang Lifted Loop WT: 2.4.4

// Common/Helper Functions 

double WT_scale_common(
    const int n
) {
    double p_1 = std::cos((2 * M_PI) / n);
    double p_2 = (3 / 8) + ((1 / 4) * p_1); 
    
    return std::pow(p_2, 2);
}

/**
 * Lifting I. 
 * The first lifting step subtracts a filtered version of the boundary vertices in Vnew from those in Vold.
 * Let v be a boundary vertex in Vold, and v1 and v2 be boundary vertices from Vnew. The relative positions of v1
 * and v2 to v are specified by the mask shown in Figure 2.16(a). Then, the updated value v` of v is given by
*/

void WT_Lifting_1(
    Eigen::Vector3d& v,
    Eigen::Vector3d v1,
    Eigen::Vector3d v2,
    Eigen::Vector3d& v_prime
) {
    v_prime = v - (0.25 * (v1 + v2));
}

/**
 * Lifting II. 
 * The second lifting step subtracts a filtered version of the boundary vertices in Vold from those in
 * Vnew. Let v be a boundary vertex in Vnew, and v1 and v2 be vertices from Vold. The relative positions of v1 and v2
 * to v are specified by the mask shown in Figure 2.16(b). Then, the updated value v` of v is given by
*/

void WT_Lifting_2(
    Eigen::Vector3d& v,
    Eigen::Vector3d v1,
    Eigen::Vector3d v2,
    Eigen::Vector3d& v_prime
) {
    v_prime = v - (0.5 * (v1 + v2));
}

// Scalar delta function for Lifting III (3)

double scalar_delta(
    const int n
) {
    double n_1 = 1 / n;
    double p_4 = 1 - ((8 / 5) * WT_scale_common(n));

    return n_1 * p_4;
}

/**
 * Lifting III. 
 * The third lifting step subtracts a filtered version of the interior vertices in Vnew from those in Vold.
 * Let v be a vertex in Vold, and v1, v2 ..., and vn be the n 1-ring neighbours of v, which are also in Vnew. The
 * relative positions of v1, v2 ..., and vn to v are specified by the mask shown in Figure 2.16(c). Then, the updated
 * value v` of v is given by
*/

void WT_Lifting_3(
    const Eigen::MatrixX3d vertices,
    const Eigen::Vector3d v,
    Eigen::Vector3d& v_prime
) {
    int len = vertices.rows();
    
    Eigen::Vector3d sum;
    sum << 0, 0, 0;

    for(int j = 0; j < len; j++) {
        sum += vertices.row(j);
    }

    v_prime = v - (scalar_delta(len) * sum);
}


/**
 * Scaling. 
 * The scaling step divides the interior vertices in Vold by a scalar. Let v be an interior vertex of valence
 * n in Vold. The updated value of v` of v is given by
*/

void WT_Scaling(
    const Eigen::Vector3d v,
    const double n,
    Eigen::Vector3d& v_prime 
) {
    v_prime = v / ((8 / 5) * WT_scale_common(n));    
}

/**
 * Lifting IV. 
 * The fourth lifting step subtracts a filtered version of the vertices in Vold from the interior vertices in
 * Vnew. Let v be an interior vertex in Vnew, and v1, v2, v3, and v4 be vertices from Vold. The relative positions of v1,
*/

void WT_Lifting_4(
    const Eigen::Vector3d v,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    const Eigen::Vector3d v3,
    const Eigen::Vector3d v4,
    Eigen::Vector3d& v_prime
) {
    const Eigen::Vector3d three_eigths = (3/8) * (v1 + v2);
    const Eigen::Vector3d one_eigth = (1/8) * (v3 + v4); 

    v_prime = v - (three_eigths + one_eigth); 
}

/** 
 * Lifting V. 
 * The fifth lifting step subtracts a filtered version of the boundary vertices in Vnew from those in Vold
 * Let v1, v2, v3, and v4 be boundary vertices in Vold, and v be a boundary vertex from Vnew. The relative positions
 * of v1, v2, v3, and v4 to v are specified by the mask shown in Figure 2.16(e). Then, the updated values of v`1, v`2, v`3,
 * and v`4 of v1, v2, v3, and v4 are given by
*/

void WT_Lifting_5(
    const Eigen::Vector3d v,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    const Eigen::Vector3d v3,
    const Eigen::Vector3d v4,
    Eigen::Vector3d& v1_prime,
    Eigen::Vector3d& v2_prime,
    Eigen::Vector3d& v3_prime,
    Eigen::Vector3d& v4_prime
) {
    const double n_14 = -0.525336;
    const double n_23 = 0.189068;

    v1_prime = v1 - (n_14 * v);
    v2_prime = v2 - (n_23 * v);
    v3_prime = v3 - (n_23 * v);
    v4_prime = v4 - (n_14 * v);
}

/**
 * Lifting VI
*/

double WT_Coefficient_Alpha(
    const double n
) {
    return (3 / 8) + WT_scale_common(n); 
}

double WT_Coefficient_Beta(
    const double n
) {
    return (8 / 5) * WT_scale_common(n); 
}

double WT_Coefficient_Gamma(
    const double n
) {
    return (1 / n) * ((5 / 8) - WT_scale_common(n)); 
}

double WT_Coefficient_Delta(
    const double n
) {
    return (1 / n) * (1 - ((8 / 5) * WT_scale_common(n))); 
}