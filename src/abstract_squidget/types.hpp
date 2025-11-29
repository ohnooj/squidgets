#ifndef TYPES_HPP
#define TYPES_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <maya/MMatrix.h>

#include <maya/MDagPath.h>
#include <maya/MMatrix.h>
#include <maya/MPointArray.h>
#include <maya/MString.h>
#include <maya/MStringArray.h>
#include <maya/MVector.h>

// S P C W L o
namespace abstract_squidget {
struct RenderMatrices {
  MMatrix localMatrix;
  MMatrix worldMatrix;
  MMatrix cameraMatrix;
  MMatrix projMatrix;
  MMatrix screenMatrix;
  //   MMatrix screenTransformMatrix;
};

using V2 = Eigen::Vector2d;
using V3 = Eigen::Vector3d;
using V4 = Eigen::Vector4d;
using VX = Eigen::VectorXd;
using VXi = Eigen::VectorXi;

using M2x2 = Eigen::Matrix2d;
using M3x3 = Eigen::Matrix3d;
using M4x4 = Eigen::Matrix4d;

using M2xX = Eigen::Matrix2Xd;
using MXx2 = Eigen::MatrixX2d;
using M3xX = Eigen::Matrix3Xd;
using MXx3 = Eigen::MatrixX3d;
using M4xX = Eigen::Matrix4Xd;
using MXxX = Eigen::MatrixXd;

struct EditFlags {
  bool deform;
  bool preSelect;
  bool translateOnly;
  bool isRigCurves;
};

struct CurveDistanceItem {
  double shape_sim;
  double matrix_deviation;
  RenderMatrices renderMats;
  RenderMatrices initRenderMats;
  MDagPath geomPath;
  MPointArray squidgetCurve;
  MPointArray shadowCurve;
  double score;

  CurveDistanceItem() { score = -1; }
  CurveDistanceItem(double shape_similar, double matrix_deviation,
                    RenderMatrices renderMats, RenderMatrices initRenderMats,
                    MDagPath geomPath, MPointArray squidgetCurve,
                    MPointArray shadowCurve)
      : shape_sim(shape_similar), matrix_deviation(matrix_deviation),
        renderMats(renderMats), initRenderMats(initRenderMats),
        geomPath(geomPath), squidgetCurve(squidgetCurve),
        shadowCurve(shadowCurve) {
    // Higher score comes first

    // Normalize the distance score
    double max_dist = 200;
    double normalized_distance = std::min(matrix_deviation / max_dist, 1.0);

    // Normalize the shape similarity
    double normalized_shape_sim = std::exp(-shape_similar / 9);

    double alpha = std::pow(1 - normalized_distance, 2);

    score = alpha * 2 + (1 - alpha) * normalized_shape_sim;

    // Old scoring
    // double alpha = 0.65;
    // score = (1 - alpha) * 1 / shape_sim + alpha * 1 / matrix_deviation;
  }
};

struct CompareCurveDistanceItem {
  bool operator()(const CurveDistanceItem &lhs,
                  const CurveDistanceItem &rhs) const {
    double lhs_score = lhs.score;
    double rhs_score = rhs.score;
    return lhs_score < rhs_score;
  }
};

struct QueryClosestObjectResult {
  MString object;
  MPointArray preShadowCurve;
  MPointArray postShadowCurve;
  RenderMatrices renderMats;
  RenderMatrices initRenderMats;
};

enum DeformType { RIGID_TRANSFORM, CURVE_DEFORM };
struct DeformObjectResult {
  DeformType deformType;
  // Rigid Transform
  MTransformationMatrix oldTransform;
  MVector offset;
  // Curve Deform
  MPointArray oldVertices;
  MPointArray postShadowCurve;
};

} // namespace abstract_squidget
#endif
