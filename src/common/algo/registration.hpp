#ifndef MATH_HPP
#define MATH_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <vector>

namespace common {
namespace algo {

namespace {
using namespace Eigen;
using namespace std;
using namespace cv;
} // namespace

struct Edge {
  int point1;
  int point2;
};

double euclidean_distance(const VectorXd &point1, const VectorXd &point2);

typedef Vector2d Point; // Define a 2D point using Eigen
tuple<double, double, double>
point_based_matching(const vector<pair<Point, Point>> &point_pairs);

vector<pair<double, int>> nearest_neighbors(const MatrixX2d &reference_points,
                                            const MatrixX2d &points);

pair<vector<Matrix3d>, MatrixX2d>
icp(const MatrixX2d &reference_points, MatrixX2d points,
    int max_iterations = 100, double distance_threshold = 0.3,
    double convergence_translation_threshold = 1e-3,
    double convergence_rotation_threshold = 1e-4,
    int point_pairs_threshold = 10, bool verbose = false);

// Point to Plane ICP ==========================================================
void random_points_on_mesh(const int n, const MatrixXd &V, const MatrixXi &E,
                           MatrixX2d &X);

void calculate_edge_normals(const MatrixX2d &VX, const MatrixX2i &EX,
                            MatrixX2d &N);

void calculate_vertex_normals(const MatrixX2d &VX, const MatrixX2i &EX,
                              const MatrixX2d &EN, MatrixX2d &VN);

void calculate_interpolated_normals(const MatrixX2d &Y, const MatrixX2i &EY,
                                    const MatrixX2d &VNY, const MatrixX2d &C,
                                    const VectorXi &I, MatrixX2d &C_normals);

VectorXd compute_edge_lengths(const MatrixXd &V, const MatrixXi &E);

VectorXd cumulative_sum(const VectorXd &input);

void point_to_plane_icp_open(const MatrixX2d &_X, const MatrixX2d &Y,
                             const int number_samples, const int max_iterations,
                             Matrix2d &R, RowVector2d &t,
                             vector<MatrixX2d> &sample_list);

MatrixX2d findClosestPoints(const MatrixX2d &src, const MatrixX2d &targ);
// Estimate affine transformation using OpenCV  ================================
void estimate_affine_partial(const MatrixX2d &_X, const MatrixX2d &Y,
                             const int number_samples, const int max_iterations,
                             Matrix2d &R, RowVector2d &t,
                             vector<MatrixX2d> &sample_list);

// Function to compute the centroid of a set of points
Eigen::Vector2d computeCentroid(const Eigen::MatrixX2d &points);

// Function to subtract the centroid from a set of points
Eigen::MatrixX2d subtractCentroid(const Eigen::MatrixX2d &points,
                                  const Eigen::Vector2d &centroid);

// Function to compute the covariance matrix
Eigen::Matrix2d computeCovarianceMatrix(const Eigen::MatrixX2d &points1,
                                        const Eigen::MatrixX2d &points2);

// Function to find the rotation matrix using SVD
Eigen::Matrix2d findRotationMatrix(const Eigen::Matrix2d &covariance);

// Function to compute the rigid transformation (rotation and translation)
void rigidTransformation(const Eigen::MatrixX2d &points1,
                         const Eigen::MatrixX2d &points2,
                         Eigen::Matrix2d &rotation,
                         Eigen::RowVector2d &translation);

} // namespace algo
} // namespace common
#endif