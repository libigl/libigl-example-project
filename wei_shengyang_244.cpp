#include <igl/PI.h>
#include <igl/upsample.h>

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
    return (1 / n) * (1 - (8 / 5) * WT_scale_common(n)); 
}

void WT_Get_A(
    const double n_0,
    const double n_1,
    const double n_2,
    const double n_3,
    Eigen::MatrixX4d& A
) {
    // Coefficients

    const double alpha_1 = WT_Coefficient_Alpha(1);
    const double alpha_2 = WT_Coefficient_Alpha(2);
    const double alpha_3 = WT_Coefficient_Alpha(3);
    const double alpha_4 = WT_Coefficient_Alpha(4);

    const double gamma_0 = WT_Coefficient_Gamma(0);
    const double gamma_1 = WT_Coefficient_Gamma(1);
    const double gamma_2 = WT_Coefficient_Gamma(2);
    const double gamma_3 = WT_Coefficient_Gamma(3);
    const double gamma_4 = WT_Coefficient_Gamma(4);

    /** START OF VARIABLE INITIALIZATION **/

    /** 
     * IDENTITY BOYOS
    */
    double a_11 = 0;
    double a_22 = 0;
    double a_33 = 0;
    double a_44 = 0;

    /** 
     * Mirror-able
     * These values can be mirrored. i.e. a_12 == a_21, a_13 == a_31, etc...
    */
    double a_12 = 0;
    double a_13 = 0;
    double a_14 = 0;
    double a_23 = 0;
    double a_24 = 0;
    double a_34 = 0;

    /** END OF VARIABLE INITIALIZATION **/

    // Hair-loss-math-time (its not hard just exhaustive)
    // first start with the IDENTITY BOYOS

    a_11 += std::pow(alpha_1, 2);
    a_11 += std::pow(gamma_2, 2);
    a_11 += std::pow(gamma_3, 2);
    a_11 += std::pow(gamma_4, 2);
    a_11 += (1 / 256) * (n_0 - 3);
    a_11 += (5 / 32) * n_0;

    a_22 += std::pow(gamma_1, 2);
    a_22 += std::pow(alpha_2, 2);
    a_22 += std::pow(gamma_3, 2);
    a_22 += std::pow(gamma_4, 2);
    a_22 += (1 / 256) * (n_1 - 3);
    a_22 += (5 / 32) * n_1;

    a_33 += std::pow(gamma_0, 2);
    a_33 += std::pow(gamma_1, 2);
    a_33 += std::pow(alpha_2, 2);
    a_33 += (1 / 256) * (n_2 - 3);
    a_33 += (5 / 32) * n_2;

    a_44 += std::pow(gamma_0, 2);
    a_44 += std::pow(gamma_1, 2);
    a_44 += std::pow(alpha_3, 2);
    a_44 += (1 / 256) * (n_3 - 3);
    a_44 += (5 / 32) * n_3;

    // lastly the mirror bois

    a_12 += alpha_1 * gamma_1;
    a_12 += gamma_2 * alpha_1;
    a_12 += std::pow(gamma_3, 2);
    a_12 += std::pow(gamma_4, 2);
    a_12 += (21 / 64);

    a_13 += alpha_1 * gamma_1;
    a_13 += std::pow(gamma_2, 2);
    a_13 += gamma_3 * alpha_3;
    a_13 += (85 / 256);

    a_14 += alpha_1 * gamma_1;
    a_14 += std::pow(gamma_2, 2);
    a_14 += gamma_4 * alpha_4;
    a_14 += (85 / 256);
    
    a_23 += std::pow(gamma_0, 2);
    a_23 += alpha_1 * gamma_1;
    a_23 += gamma_2 * alpha_2;
    a_23 += (85 / 256);

    a_24 += std::pow(gamma_0, 2);
    a_24 += alpha_1 * gamma_1;
    a_24 += gamma_3 * alpha_3;
    a_24 += (85 / 256);

    a_34 += std::pow(gamma_0, 2);
    a_34 += std::pow(gamma_1, 2);
    a_34 += (1 / 64);

    A << 
        a_11, a_12, a_13, a_14,
        a_12, a_22, a_23, a_24,
        a_13, a_23, a_33, a_34,
        a_14, a_24, a_34, a_44;
}

void WT_GET_B(
    Eigen::Vector4d& B
) {
    const double alpha_0 = WT_Coefficient_Alpha(0);
    const double alpha_1 = WT_Coefficient_Alpha(1);

    const double gamma_0 = WT_Coefficient_Gamma(0);
    const double gamma_1 = WT_Coefficient_Gamma(1);
    
    const double delta_0 = WT_Coefficient_Delta(0);
    const double delta_1 = WT_Coefficient_Delta(1);

    double b_1 = 0;
    double b_2 = 0;
    double b_3 = 0;

    b_1 += alpha_0 * delta_0;
    b_1 += gamma_1 * delta_1;
    b_1 += (3 / 8);

    b_2 += gamma_0 * delta_0;
    b_2 += alpha_1 * delta_1;
    b_2 += (3 / 8);

    b_3 += gamma_0 * delta_0;
    b_3 += gamma_1 * delta_1;
    b_3 += (1 / 8);

    B << 
        b_1,
        b_2,
        b_3,
        b_3;
}

void WT_Solve_Weights(
    const double n_0,
    const double n_1,
    const double n_2,
    const double n_3,
    Eigen::Vector4d& W
) {
    Eigen::MatrixX4d A;
    WT_Get_A(n_0, n_1, n_2, n_3, A);

    Eigen::Vector4d B;
    WT_GET_B(B);

    W = A.inverse() * B;
}   

