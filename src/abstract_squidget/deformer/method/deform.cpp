#include "abstract_squidget/deformer/method/deform.hpp"

namespace abstract_squidget {
namespace deformer {
namespace method {
void findInterval(MPointArray meshVertices, MPointArray shadowCurve,
                  int &start_interval, int &end_interval, int &increment) {
  /**
   * @brief Find the interval of mesh vertices that are closest to the shadow
   * curve.
   */
  Eigen::MatrixX2d X(meshVertices.length(), 2);
#pragma omp parallel for
  for (int i = 0; i < meshVertices.length(); ++i) {
    X(i, 0) = meshVertices[i].x;
    X(i, 1) = meshVertices[i].y;
  }

  Eigen::MatrixX2i EX(X.rows(), 2);
#pragma omp parallel for
  for (int i = 0; i < EX.rows(); ++i) {
    EX(i, 0) = i;
    EX(i, 1) = (i + 1) % X.rows();
  }

  igl::AABB<Eigen::MatrixX2d, 2> xtree;
  xtree.init(X, EX);

  Eigen::MatrixX2d Y(shadowCurve.length(), 2);
  for (size_t i = 0; i < shadowCurve.length(); ++i) {
    Y(i, 0) = shadowCurve[i].x;
    Y(i, 1) = shadowCurve[i].y;
  }

  Eigen::VectorXd D;
  Eigen::VectorXi I;
  Eigen::MatrixX2d C;
  xtree.squared_distance(X, EX, Y, D, I, C);

  start_interval = EX(I(0), 0);
  end_interval = EX(I(I.rows() - 1), 0);

  // Determine the direction the interval is heading
  // If more pairs of I(i) I(i+1) are ascending, then increment is positive
  int increase_count = 0;
  int decrease_count = 0;
  for (int i = 0; i < I.rows() - 1; ++i) {
    int a = EX(I(i), 0);
    int b = EX(I(i + 1), 0);
    if (a == 0 && b == int(meshVertices.length() - 1)) {
      decrease_count += 1;
    } else if (a == int(meshVertices.length() - 1) && b == 0) {
      increase_count += 1;
    } else if (a < b) {
      increase_count += 1;
    } else {
      decrease_count += 1;
    }
  }

  increment = (increase_count > decrease_count) ? 1 : -1;
}

Eigen::VectorXd parameterize_curve(const Eigen::MatrixX2d &points) {
  int n = points.rows();
  Eigen::VectorXd cum_dist(n);
  cum_dist(0) = 0;

  // Calculate distances between consecutive points
  for (int i = 1; i < n; ++i) {
    cum_dist(i) = cum_dist(i - 1) + (points.row(i) - points.row(i - 1)).norm();
  }

  // Normalize distances to range from 0 to 1
  double max_distance = cum_dist(n - 1);
  Eigen::VectorXd normalized_cum_dist = cum_dist / max_distance;
  return normalized_cum_dist;
}

Eigen::MatrixX2d convertMPointArrayToMatrixX2d(const MPointArray &mPointArray) {
  // Initialize an Eigen matrix with the same number of rows as mPointArray and
  // 2 columns for x, y
  Eigen::MatrixX2d matrix(mPointArray.length(), 2);

  // Populate the matrix
#pragma omp parallel for
  for (int i = 0; i < mPointArray.length(); ++i) {
    matrix(i, 0) = mPointArray[i].x; // X coordinate
    matrix(i, 1) = mPointArray[i].y; // Y coordinate
  }

  return matrix;
}

void convertPointsToScreenSpace(const MPointArray &objSpaceV,
                                const RenderMatrices &renderMats,
                                MPointArray &screenSpaceV) {
  screenSpaceV.clear();
  screenSpaceV.setLength(objSpaceV.length());
#pragma omp parallel for
  for (int i = 0; i < objSpaceV.length(); ++i) {
    MPoint screenp;
    convertPointToScreenSpace(objSpaceV[i], renderMats, screenp);
    screenSpaceV.set(screenp, i);
  }
}

void convertPointToScreenSpace(const MPoint &objSpaceV,
                               const RenderMatrices &renderMats,
                               MPoint &screenSpaceV) {
  MPoint screenp = renderMats.screenMatrix * renderMats.projMatrix *
                   renderMats.cameraMatrix * renderMats.worldMatrix *
                   renderMats.localMatrix * objSpaceV;
  screenp.cartesianize();
  screenSpaceV = screenp;
}

Eigen::VectorXd findPointOnCurve(const Eigen::MatrixXd &Y,
                                 const Eigen::VectorXd &Y_t, double t) {
  int n = Y_t.rows();
  if (n == 0) {
    throw std::runtime_error("The input Y_t vector is empty.");
  }

  // Binary search to find the greatest element less than t
  int low = 0, high = n;
  while (low < high) {
    int mid = low + (high - low) / 2;
    if (Y_t(mid) < t) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }

  // `low` is now the index of the first element greater than or equal to `t`
  // We need to interpolate between `low - 1` (if it exists) and `low`
  int index = low - 1; // This is the greatest element less than `t`

  // Handle edge cases
  if (index < 0) {
    // t is less than the smallest element in Y_t
    return Y.row(0);
  } else if (index >= n - 1) {
    // t is greater than or equal to the last element in Y_t
    return Y.row(n - 1);
  } else {
    // Perform linear interpolation between Y(index) and Y(index + 1)
    double t0 = Y_t(index);
    double t1 = Y_t(index + 1);
    double factor = (t - t0) / (t1 - t0);
    return Y.row(index) + factor * (Y.row(index + 1) - Y.row(index));
  }
}
} // namespace method
} // namespace deformer
} // namespace abstract_squidget