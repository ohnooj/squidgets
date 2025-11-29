#include "abstract_squidget/deformer/optimization.hpp"

#include "abstract_squidget/deformer/deformerCurves.hpp"
#include "abstract_squidget/deformer/renderMatrices.hpp"
#include "common/algo/registration.hpp"
#include "common/util/maya.hpp"

namespace abstract_squidget {
namespace deformer {
RenderMatrices optimizeNewMatrices(MString object, MPointArray strokePoints2D,
                                   MPointArray shadowCurve, EditFlags flags) {
  MDagPath geomPath = common::util::getDagPathFromName(object);
  RenderMatrices renderMats = getRenderMatrices(geomPath.node());

  renderMats =
      optimizeNewMatrices(shadowCurve, strokePoints2D, renderMats, flags);
  return renderMats;
}

RenderMatrices optimizeNewMatrices(MPointArray squidgetObjCVs,
                                   MPointArray curvePoints2D,
                                   RenderMatrices renderMats, EditFlags flags) {
  cout << "optimization::optimizeNewMatrices()" << endl;
  if (squidgetObjCVs.length() < 2) {
    return renderMats;
  }
  // A S P C W L O = S P C W B L O
  // A X L O = X B L O
  // A X = X B
  // X^-1 A X = B
  M4x4 A = leastSquares(squidgetObjCVs, curvePoints2D, renderMats, flags);
  M4x4 X = common::util::mMatrixToMatrix(renderMats.screenMatrix) *
           common::util::mMatrixToMatrix(renderMats.projMatrix) *
           common::util::mMatrixToMatrix(renderMats.cameraMatrix) *
           common::util::mMatrixToMatrix(renderMats.worldMatrix);
  M4x4 B = X.inverse() * A * X;

  MMatrix newLocalMatrix =
      common::util::matrixToMMatrix(B) * renderMats.localMatrix;

  RenderMatrices ret;
  ret.localMatrix = newLocalMatrix;
  ret.worldMatrix = renderMats.worldMatrix;
  ret.cameraMatrix = renderMats.cameraMatrix;
  ret.projMatrix = renderMats.projMatrix;
  ret.screenMatrix = renderMats.screenMatrix;
  return ret;
}

// =============================================================================
//                            Least Squares
// =============================================================================
Eigen::Matrix4d leastSquares(MPointArray initialPoints,
                             MPointArray targetPoints,
                             RenderMatrices renderMats, EditFlags flags) {
  /**
   * @brief Perform Least Squares by projecting <initialPoints> to screen
   * space and finding the transformation matrix that maps the project
   * <initialPoints> to <targetPoints>.
   */

  // Project the initial points to screen space
  M4xX o = common::util::pointArrayToMatrix(initialPoints);
  M4xX objScreenPts = getScreenPoints(o, renderMats);

  // Convert all 3D points to 2D homogeneous points
  M3xX objScreen2d = objScreenPts.block(0, 0, 3, objScreenPts.cols());
  objScreen2d.row(2) = VX::Ones(objScreen2d.cols());

  // Project the target points to screen space
  M4xX v = common::util::pointArrayToMatrix(targetPoints);
  M3xX v2d = v.block(0, 0, 3, v.cols());
  v2d.row(2) = VX::Ones(v2d.cols());

  M4x4 screenA4 = M4x4::Identity();
  if (v2d.cols() < 2 || objScreen2d.cols() < 2) {
  // if (v2d.cols() < 3 || objScreen2d.cols() < 3) {
    cout << "Not enough points to optimize" << endl;
    return screenA4;
  }

  if (flags.translateOnly) {
    screenA4 = optimizeTranslation(v2d, objScreen2d, flags);
  } else {
    Eigen::MatrixX2d query_curve = v.block(0, 0, 2, v.cols()).transpose();
    Eigen::MatrixX2d shape_curve =
        objScreen2d.block(0, 0, 2, objScreen2d.cols()).transpose();

    Eigen::Matrix2d currR;
    Eigen::RowVector2d currT;
    std::vector<Eigen::MatrixX2d> sample_list;

    if (flags.preSelect) {
      // Use point-to-point only if we have a preselect curve.
      // common::algo::estimate_affine_partial(query_curve, shape_curve, 60, 4,
      //                                       currR, currT, sample_list);
      common::algo::point_to_plane_icp_open(query_curve, shape_curve, 60, 3,
                                            currR, currT, sample_list);
    } else {
      common::algo::point_to_plane_icp_open(query_curve, shape_curve, 60, 5,
                                            currR, currT, sample_list);
    }

    screenA4 = M4x4::Identity();
    screenA4.block(0, 0, 2, 2) = currR;
    screenA4.block(3, 0, 1, 2) = currT;
    screenA4.transposeInPlace();
    M4x4 screenA4inv = screenA4.inverse();
    screenA4 = screenA4inv;
  }

  return screenA4;
}

// =============================================================================
//                            Optimize Translation
// =============================================================================

// Eigen::Matrix4d optimizeTranslation(M3xX v2d, M3xX objScreen2d,
//                                     EditFlags flags) {
//   // Optimize Translation between v2d and objScreen2d.shadow
//   // cout << "optimization::optimizeTranslation()" << endl;
//   Eigen::MatrixX2d _objScreenCoords =
//       objScreen2d.block(0, 0, 2, objScreen2d.cols()).transpose();
//   Eigen::MatrixX2d strokeScreenCoords =
//       v2d.block(0, 0, 2, v2d.cols()).transpose();

//   Eigen::MatrixX2d objScreenCoords = _objScreenCoords;
//   if (!flags.preSelect) {
//     // Use point-to-point only if we have a preselect curve.
//     objScreenCoords =
//         algo::findClosestPoints(strokeScreenCoords, _objScreenCoords);
//   }

//   // Find Centroids of both sets of points
//   Eigen::Vector2d objScreenCentroid = algo::computeCentroid(objScreenCoords);
//   Eigen::Vector2d strokeScreenCentroid =
//       algo::computeCentroid(strokeScreenCoords);

//   // Create Matrix using only translation of centroids
//   Eigen::Vector2d centroidDiff = strokeScreenCentroid - objScreenCentroid;
//   M4x4 translation = M4x4::Identity();
//   translation.block(0, 3, 2, 1) = centroidDiff;

//   return translation;
// }

Eigen::Matrix4d optimizeTranslation(M3xX v2d, M3xX objScreen2d,
                                    EditFlags flags) {
  // Optimize Translation between v2d and objScreen2d.shadow
  // cout << "optimization::optimizeTranslation()" << endl;
  Eigen::MatrixX2d _objScreenCoords =
      objScreen2d.block(0, 0, 2, objScreen2d.cols()).transpose();
  Eigen::MatrixX2d strokeScreenCoords =
      v2d.block(0, 0, 2, v2d.cols()).transpose();

  Eigen::MatrixX2d objScreenCoords = _objScreenCoords;

  Eigen::Vector2d T = Eigen::Vector2d::Zero();
  for (int it = 0; it < 1; it++) {
    if (!flags.preSelect) {
      // Use point-to-point only if we have a preselect curve.
      objScreenCoords = common::algo::findClosestPoints(
          strokeScreenCoords, _objScreenCoords.rowwise() + T.transpose());
    }

    // Find Centroids of both sets of points
    Eigen::Vector2d objScreenCentroid =
        common::algo::computeCentroid(objScreenCoords);
    Eigen::Vector2d strokeScreenCentroid =
        common::algo::computeCentroid(strokeScreenCoords);

    // Create Matrix using only translation of centroids
    Eigen::Vector2d centroidDiff = strokeScreenCentroid - objScreenCentroid;
    T += centroidDiff;

    // Apply Translation to objScreenCoords
    objScreenCoords = _objScreenCoords.rowwise() + T.transpose();
  }

  M4x4 translation = M4x4::Identity();
  translation.block(0, 3, 2, 1) = T;

  return translation;
}

} // namespace deformer
} // namespace abstract_squidget