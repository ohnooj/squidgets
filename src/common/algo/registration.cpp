#include "common/algo/registration.hpp"

#include <iostream>
#include <random>

#include <Eigen/LU>
#include <Eigen/SVD>
#include <igl/AABB.h>
#include <igl/rigid_alignment.h>
#include <opencv2/surface_matching/icp.hpp>

namespace common {
namespace algo {
// ############################################################################
//                          Point to Point ICP
// ############################################################################
double euclidean_distance(const VectorXd &point1, const VectorXd &point2) {
  VectorXd diff = point1 - point2;
  return diff.norm();
}

tuple<double, double, double>
point_based_matching(const vector<pair<Point, Point>> &point_pairs) {
  int n = point_pairs.size();
  if (n == 0) {
    return make_tuple(NAN, NAN, NAN);
  }

  Point sum1 = Point::Zero(), sum2 = Point::Zero();
  for (const auto &pair : point_pairs) {
    sum1 += pair.first;
    sum2 += pair.second;
  }

  Point mean1 = sum1 / n;
  Point mean2 = sum2 / n;

  double s_x_xp = 0, s_y_yp = 0, s_x_yp = 0, s_y_xp = 0;
  for (const auto &pair : point_pairs) {
    Point d1 = pair.first - mean1;
    Point d2 = pair.second - mean2;

    s_x_xp += d1.x() * d2.x();
    s_y_yp += d1.y() * d2.y();
    s_x_yp += d1.x() * d2.y();
    s_y_xp += d1.y() * d2.x();
  }

  double rot_angle = atan2(s_x_yp - s_y_xp, s_x_xp + s_y_yp);
  double translation_x =
      mean2.x() - (mean1.x() * cos(rot_angle) - mean1.y() * sin(rot_angle));
  double translation_y =
      mean2.y() - (mean1.x() * sin(rot_angle) + mean1.y() * cos(rot_angle));

  return make_tuple(rot_angle, translation_x, translation_y);
}

vector<pair<double, int>> nearest_neighbors(const MatrixX2d &reference_points,
                                            const MatrixX2d &points) {
  vector<pair<double, int>> indices(points.rows());
  for (int i = 0; i < points.rows(); ++i) {
    double min_distance = numeric_limits<double>::max();
    int min_index = -1;
    for (int j = 0; j < reference_points.rows(); ++j) {
      Point point{points(i, 0), points(i, 1)};
      Point ref_point{reference_points(j, 0), reference_points(j, 1)};
      double distance = euclidean_distance(point, ref_point);
      if (distance < min_distance) {
        min_distance = distance;
        min_index = j;
      }
    }
    indices[i] = make_pair(min_distance, min_index);
  }
  return indices;
}

pair<vector<Matrix3d>, MatrixX2d> icp(const MatrixX2d &reference_points,
                                      MatrixX2d points, int max_iterations,
                                      double distance_threshold,
                                      double convergence_translation_threshold,
                                      double convergence_rotation_threshold,
                                      int point_pairs_threshold, bool verbose) {

  vector<Matrix3d> transformation_history;
  // Implement the ICP loop with Eigen and nearest neighbor search

  for (int iter_num = 0; iter_num < max_iterations; ++iter_num) {
    if (verbose) {
      cout << "------ Iteration " << iter_num << endl;
    }

    vector<pair<double, int>> closest_point_distances = nearest_neighbors(
        reference_points,
        points); // for each point, get index of closest reference point

    vector<pair<Point, Point>> closest_point_pairs;
    for (int i = 0; i < points.rows(); ++i) {
      if (closest_point_distances[i].first < distance_threshold) {
        int cp_index = closest_point_distances[i].second;
        auto first_point = Point(points(i, 0), points(i, 1));
        auto second_point =
            Point(reference_points(cp_index, 0), reference_points(cp_index, 1));
        closest_point_pairs.emplace_back(first_point, second_point);
      }
    }

    if (verbose) {
      cout << "Found " << closest_point_pairs.size() << " point pairs." << endl;
    }
    if ((int)closest_point_pairs.size() < point_pairs_threshold) {
      cout << "Not enough point pairs found. Stopping." << endl;
      break;
    }

    tuple<double, double, double> result =
        point_based_matching(closest_point_pairs);
    double closest_rot_angle = get<0>(result);
    double closest_translation_x = get<1>(result);
    double closest_translation_y = get<2>(result);

    if (closest_rot_angle != NAN) {
      if (verbose) {
        cout << "Rotation: " << closest_rot_angle << endl;
        cout << "Translation: " << closest_translation_x << ", "
             << closest_translation_y << endl;
      }
    }
    if (closest_rot_angle == NAN || closest_translation_x == NAN ||
        closest_translation_y == NAN) {
      if (verbose) {
        cout << "NAN result. Stopping." << endl;
      }
      break;
    }

    // Transform points using the calculated rotateion and translation
    double c = cos(closest_rot_angle);
    double s = sin(closest_rot_angle);

    // Create the rotation matrix
    Matrix2d rot;
    rot << c, -s, s, c;

    // Apply the rotation
    MatrixX2d rotated_points = points * rot.transpose();

    // Apply the translation
    for (int i = 0; i < rotated_points.rows(); ++i) {
      rotated_points(i, 0) += closest_translation_x;
      rotated_points(i, 1) += closest_translation_y;
    }

    points = rotated_points;

    // update transformation history
    Matrix3d transformation;
    transformation << c, -s, closest_translation_x, s, c, closest_translation_y,
        0, 0, 1;
    transformation_history.push_back(transformation);

    // Check for convergence
    if (abs(closest_translation_x) < convergence_translation_threshold &&
        abs(closest_translation_y) < convergence_translation_threshold &&
        abs(closest_rot_angle) < convergence_rotation_threshold) {
      if (verbose) {
        cout << "Converged." << endl;
      }
      break;
    }
  }

  // Return the transformation history and the final aligned points
  return {transformation_history, points};
}

// ############################################################################
//                          Point to Plane ICP
// ############################################################################
/**
 * @brief Compute the edge normals for the given vertices and edges.
 *
 * @param VX
 * @param EX
 * @param N output normals
 */
void calculate_edge_normals(const MatrixX2d &VX, const MatrixX2i &EX,
                            MatrixX2d &N) {
  for (int i = 0; i < EX.rows(); ++i) {
    RowVector2d edge = VX.row(EX(i, 1)) - VX.row(EX(i, 0));
    RowVector2d normal = RowVector2d(edge.y(), -edge.x());
    normal.normalize();
    N.row(i) = normal;
  }
}

/**
 * @brief Compute the normals for each point in the mesh. End points have
 * normals of -edge.
 *
 * @param VX
 * @param EX
 * @param EN
 * @param VN
 */
void calculate_vertex_normals(const MatrixX2d &VX, const MatrixX2i &EX,
                              const MatrixX2d &EN, MatrixX2d &VN) {
  VN.resize(VX.rows(), 2);
  VN.setZero();

  // Create a vector to count the number of edges connected to each vertex
  vector<int> counts(VX.rows(), 0);

  // Sum the normals of the edges connected to each vertex
  for (int i = 0; i < EX.rows(); ++i) {
    int v1 = EX(i, 0);
    int v2 = EX(i, 1);
    VN.row(v1) += EN.row(i);
    VN.row(v2) += EN.row(i);
    counts[v1]++;
    counts[v2]++;
  }

  // Average the normals and handle vertices connected to only one edge
  for (int i = 0; i < VX.rows(); ++i) {
    if (counts[i] == 1) {
      for (int j = 0; j < EX.rows(); ++j) {
        if (EX(j, 0) == i) {
          int other = EX(j, 1);
          Vector2d normal = VX.row(i) - VX.row(other);
          normal.normalize();
          VN.row(i) = normal;
          break;
        } else if (EX(j, 1) == i) {
          int other = EX(j, 0);
          Vector2d normal = VX.row(i) - VX.row(other);
          normal.normalize();
          VN.row(i) = normal;
          break;
        }
      }
    } else if (counts[i] > 1) {
      VN.row(i) /= counts[i];
      VN.row(i).normalize();
    }
  }
}

void calculate_interpolated_normals(const MatrixX2d &Y, const MatrixX2i &EY,
                                    const MatrixX2d &VNY, const MatrixX2d &C,
                                    const VectorXi &I, MatrixX2d &C_normals) {
  C_normals.resize(C.rows(), 2);

  for (int i = 0; i < C.rows(); ++i) {
    int edge_index = I(i);
    int v1 = EY(edge_index, 0);
    int v2 = EY(edge_index, 1);

    Vector2d edge = Y.row(v2) - Y.row(v1);
    double edge_length = edge.norm();
    double t = (C.row(i) - Y.row(v1)).dot(edge) / (edge_length * edge_length);

    C_normals.row(i) = (1 - t) * VNY.row(v1) + t * VNY.row(v2);
    C_normals.row(i).normalize();
  }
}

/**
 * @brief Sample random points on the 2D mesh
 *
 * @param n number of points
 * @param V V x 2 matrix of vertices
 * @param E E x 2 matrix of edges
 * @param X n x 2 matrix of sampled points
 */
void random_points_on_mesh(const int n, const MatrixXd &V, const MatrixXi &E,
                           MatrixX2d &X) {
  // Random number generator
  random_device rd;
  mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<> dis(0.0, 1.0);

  // Area of each edge
  VectorXd A = compute_edge_lengths(V, E);

  // Cumulative sum of relative areas
  // Should be between [0, 1]
  VectorXd A_cum = cumulative_sum(A);
  A_cum /= A_cum[A_cum.size() - 1];

  X.resize(n, 2);
  for (int i = 0; i < n; ++i) {
    // Find random triangle index using binary search
    double a = dis(gen);
    int tri_index = A_cum.size() / 2;
    int left = 0;
    int right = A_cum.size() - 1;
    while (left < right) {
      if (a < A_cum[tri_index]) {
        right = tri_index;
      } else {
        left = tri_index + 1;
      }
      tri_index = (left + right) / 2;
    }

    // Find random point in triangle
    double t = dis(gen);
    RowVectorXi e = E.row(tri_index);
    RowVectorXd u = V.row(e[0]);
    RowVectorXd v = V.row(e[1]);
    RowVectorXd p = u + t * (v - u);

    X.row(i) = p;
  }
}

/**
 * @brief Computes the edge lengths for the given vertices and edges.
 *
 * @param V
 * @param E
 * @return VectorXd length of each edge
 */
VectorXd compute_edge_lengths(const MatrixXd &V, const MatrixXi &E) {
  MatrixXd S(E.rows(), V.cols());
  MatrixXd T(E.rows(), V.cols());
  for (int i = 0; i < E.rows(); ++i) {
    S.row(i) = V.row(E(i, 0));
    T.row(i) = V.row(E(i, 1));
  }

  MatrixXd differences = T - S;
  MatrixXd squared_differences = differences.array().square();
  VectorXd distances = squared_differences.rowwise().sum().array().sqrt();
  return distances;
}

VectorXd cumulative_sum(const VectorXd &input) {
  VectorXd result(input.size());
  if (input.size() == 0)
    return result;

  result(0) = input(0);
  for (int i = 1; i < input.size(); ++i) {
    result(i) = result(i - 1) + input(i);
  }
  return result;
}

MatrixX2d findClosestPoints(const MatrixX2d &src, const MatrixX2d &targ) {
  // cout << "registration::findClosestPoints()" << endl;

  int skip = 10;

  // int n = targ.rows() / skip;
  // MatrixX2d smallTarg(n, 2);
  // smallTarg = targ(Eigen::seqN(0, n, skip), Eigen::all);

  // MatrixX2i targEdges(n, 2);
  // targEdges.col(0) = VectorXi::LinSpaced(n, 0, n - 1);
  // targEdges.col(1) = (targEdges.col(0).array() + 1).eval();
  // targEdges(n - 1, 1) = 0; // Dunno how to mod in Eigen
  // // targEdges = targEdges * skip;
  // times.push_back(std::chrono::high_resolution_clock::now());

  // igl::AABB<MatrixX2d, 2> targTree;
  // targTree.init(smallTarg, targEdges);

  // times.push_back(std::chrono::high_resolution_clock::now());

  // VectorXd D;
  // VectorXi I;
  // targTree.squared_distance(smallTarg, targEdges, src, D, I, closestPoints);
  // times.push_back(std::chrono::high_resolution_clock::now());

  // This is literally faster than setting up the AABB tree
  MatrixX2d closestPoints;
  closestPoints = MatrixX2d(src.rows(), 2);
  for (int i = 0; i < src.rows(); ++i) {
    int closestIndex = 0;
    for (int j = 0; j < targ.rows(); j = j + skip) {
      if (euclidean_distance(src.row(i), targ.row(j)) <
          euclidean_distance(src.row(i), targ.row(closestIndex))) {
        closestIndex = j;
      }
    }
    closestPoints.row(i) = targ.row(closestIndex);
  }

  return closestPoints;
}

// ############################################################################
//                          OpenCV estimateAffinePartial2D
// ############################################################################
void estimate_affine_partial(const MatrixX2d &_X, const MatrixX2d &Y,
                             const int number_samples, const int max_iterations,
                             Matrix2d &R, RowVector2d &t,
                             vector<MatrixX2d> &sample_list) {

  // std::cout << "esimate_affine_partial()" << std::endl;
  auto X(_X);
  MatrixX2i EX(X.rows() - 1, 2);
  for (int i = 0; i < EX.rows(); ++i) {
    EX(i, 0) = i;
    EX(i, 1) = (i + 1) % X.rows();
  }

  // Create edge matrix from Y
  // MatrixX2i EY(Y.rows(), 2);
  MatrixX2i EY(Y.rows() - 1, 2);
  for (int i = 0; i < EY.rows(); ++i) {
    EY(i, 0) = i;
    EY(i, 1) = (i + 1) % Y.rows();
  }

  // Build AABB tree
  igl::AABB<MatrixX2d, 2> xtree;
  xtree.init(X, EX);
  igl::AABB<MatrixX2d, 2> ytree;
  ytree.init(Y, EY);

  R = Matrix2d::Identity();
  t = Vector2d::Zero();

  MatrixX2d samples;
  VectorXd D;
  VectorXi I;
  MatrixX2d C;

  Eigen::Vector2d X_cent = computeCentroid(X);
  Eigen::Vector2d Y_cent = computeCentroid(Y);

  // move X to Y
  Eigen::Vector2d diff = Y_cent - X_cent;

  t = t + diff.transpose();
  X = (_X * R).rowwise() + t;


  for (int it = 0; it < max_iterations; ++it) {
    // Sample points on query
    random_points_on_mesh(number_samples, X, EX, samples);
    sample_list.push_back(samples);

    // Project closest points on model
    ytree.squared_distance(Y, EY, samples, D, I, C);

    Matrix2d Rup = Matrix2d::Zero();
    RowVector2d tup = RowVector2d::Zero();
    rigidTransformation(samples, C, Rup, tup);

    Rup.transposeInPlace();

    R = (R * Rup).eval();
    t = (t * Rup + tup).eval();

    // Transform X
    X = (_X * R).rowwise() + t;
  }
}

Eigen::Vector2d computeCentroid(const Eigen::MatrixX2d &points) {
  Eigen::Vector2d centroid = points.colwise().mean();
  return centroid;
}

Eigen::MatrixX2d subtractCentroid(const Eigen::MatrixX2d &points,
                                  const Eigen::Vector2d &centroid) {
  Eigen::MatrixX2d centeredPoints = points.rowwise() - centroid.transpose();
  return centeredPoints;
}

Eigen::Matrix2d computeCovarianceMatrix(const Eigen::MatrixX2d &points1,
                                        const Eigen::MatrixX2d &points2) {
  return points1.transpose() * points2;
}

Eigen::Matrix2d findRotationMatrix(const Eigen::Matrix2d &covariance) {
  Eigen::JacobiSVD<Eigen::Matrix2d> svd(covariance, Eigen::ComputeFullU |
                                                        Eigen::ComputeFullV);
  Eigen::Matrix2d rotation = svd.matrixV() * svd.matrixU().transpose();
  return rotation;
}

void rigidTransformation(const Eigen::MatrixX2d &points1,
                         const Eigen::MatrixX2d &points2,
                         Eigen::Matrix2d &rotation,
                         Eigen::RowVector2d &translation) {
  Eigen::Vector2d centroid1 = computeCentroid(points1);
  Eigen::Vector2d centroid2 = computeCentroid(points2);

  Eigen::MatrixX2d centeredPoints1 = subtractCentroid(points1, centroid1);
  Eigen::MatrixX2d centeredPoints2 = subtractCentroid(points2, centroid2);

  Eigen::Matrix2d covariance =
      computeCovarianceMatrix(centeredPoints1, centeredPoints2);

  rotation = findRotationMatrix(covariance);

  translation = (centroid2 - rotation * centroid1).transpose();
}

// ############################################################################
//                          Point to Plane ICP OpenCV
// ############################################################################
void point_to_plane_icp_open(const MatrixX2d &_X, const MatrixX2d &_Y,
                             const int number_samples, const int max_iterations,
                             Matrix2d &R, RowVector2d &t,
                             vector<MatrixX2d> &sample_list) {
  std::cout << "point_to_plane_icp_open()" << std::endl;

  auto X(_X);
  MatrixX2i EX(X.rows() - 1, 2);
  for (int i = 0; i < EX.rows(); ++i) {
    EX(i, 0) = i;
    EX(i, 1) = (i + 1) % X.rows();
  }

  int skip = (int) (_Y.rows() * 0.01);
  if (skip < 1) {
    skip = 1;
  }
  int yn = _Y.rows() / skip;
  MatrixX2d Y(yn, 2);
  Y = _Y(Eigen::seqN(0, yn, skip), Eigen::all);

  // Create edge matrix from Y
  MatrixX2i EY(yn - 1, 2);
  for (int i = 0; i < EY.rows(); ++i) {
    EY(i, 0) = i;
    EY(i, 1) = i + 1;
  }

  // Calculate normals for each point in Y
  MatrixX2d NX = MatrixX2d::Zero(EX.rows(), 2);
  calculate_edge_normals(X, EX, NX);
  MatrixX2d NVX = MatrixX2d::Zero(X.rows(), 2);
  calculate_vertex_normals(X, EX, NX, NVX);

  // Calculate normals for each point in Y
  MatrixX2d NY = MatrixX2d::Zero(EY.rows(), 2);
  calculate_edge_normals(Y, EY, NY);
  MatrixX2d NVY = MatrixX2d::Zero(Y.rows(), 2);
  calculate_vertex_normals(Y, EY, NY, NVY);

  // Build AABB tree
  igl::AABB<MatrixX2d, 2> xtree;
  xtree.init(X, EX);
  igl::AABB<MatrixX2d, 2> ytree;
  ytree.init(Y, EY);

  R = Matrix2d::Identity();
  t = Vector2d::Zero();

  MatrixX2d samples;
  VectorXd D;
  VectorXi I;
  MatrixX2d C;

  int initial_t = -2;
  for (int i = initial_t; i < 0; ++i) {

    random_points_on_mesh(number_samples, X, EX, samples);
    ytree.squared_distance(Y, EY, samples, D, I, C);

    // Calculate difference between C and samples
    // This is done because point_to_plane_rigid_matching with far points
    // doesn't work well.
    MatrixX2d diff = C - samples;
    RowVector2d mean_diff = diff.colwise().mean();
    t = t + mean_diff;

    X = (_X * R).rowwise() + t;
  }

  for (int it = 0; it < max_iterations; ++it) {
    random_points_on_mesh(number_samples, X, EX, samples);
    sample_list.push_back(samples);

    xtree.squared_distance(X, EX, samples, D, I, C);
    MatrixX2d curr_NX = MatrixX2d::Zero(C.rows(), 2);
    calculate_interpolated_normals(X, EX, NVX, C, I, curr_NX);
    Mat xPC(samples.rows(), 6, CV_32F);
    for (int i = 0; i < samples.rows(); ++i) {
      xPC.at<float>(i, 0) = samples(i, 0);
      xPC.at<float>(i, 1) = samples(i, 1);
      xPC.at<float>(i, 2) = 0;
      xPC.at<float>(i, 3) = curr_NX(i, 0);
      xPC.at<float>(i, 4) = curr_NX(i, 1);
      xPC.at<float>(i, 5) = 0;
    }

    ytree.squared_distance(Y, EY, samples, D, I, C);
    MatrixX2d curr_NY = MatrixX2d::Zero(C.rows(), 2);
    calculate_interpolated_normals(Y, EY, NVY, C, I, curr_NY);
    Mat yPC(samples.rows(), 6, CV_32F);
    for (int i = 0; i < samples.rows(); ++i) {
      yPC.at<float>(i, 0) = C(i, 0);
      yPC.at<float>(i, 1) = C(i, 1);
      yPC.at<float>(i, 2) = 0;
      yPC.at<float>(i, 3) = curr_NY(i, 0);
      yPC.at<float>(i, 4) = curr_NY(i, 1);
      yPC.at<float>(i, 5) = 0;
    }

    double residual;
    Matx44d pose;
    cv::ppf_match_3d::ICP ic(10);
    int result = ic.registerModelToScene(xPC, yPC, residual, pose);

    if (residual > 1) {
      cout << "Residual too high. Stopping: " << residual << endl;
      break;
    }

    Matrix2d Rup = Matrix2d::Zero();
    RowVector2d tup = RowVector2d::Zero();

    Rup(0, 0) = pose(0, 0);
    Rup(0, 1) = pose(1, 0);
    Rup(1, 0) = pose(0, 1);
    Rup(1, 1) = pose(1, 1);

    tup(0) = pose(0, 3);
    tup(1) = pose(1, 3);

    R = (R * Rup).eval();
    t = (t * Rup + tup).eval();

    // Transform X
    X = (_X * R).rowwise() + t;
  }
}
} // namespace algo
} // namespace common
