#ifndef BOOKMARK_SQUIDGET_TOOL_HPP
#define BOOKMARK_SQUIDGET_TOOL_HPP

#include <queue>
#include <random>
#include <string>
#include <vector>

#include <maya/MArgDatabase.h>
#include <maya/MArgList.h>
#include <maya/MDagPath.h>
#include <maya/MGlobal.h>
#include <maya/MPointArray.h>
#include <maya/MSelectionList.h>
#include <maya/MStringArray.h>
#include <maya/MSyntax.h>

#include <maya/MFnMesh.h>
#include <maya/MItDag.h>

#include <maya/MEulerRotation.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnPfxGeometry.h>
#include <maya/MFnTransform.h>

#include <maya/MItGeometry.h>
#include <maya/MQuaternion.h>

#include <maya/MPxToolCommand.h>

#include "bookmark_squidget/canvas/mayaAttributeEditor.hpp"
#include "common/algo/registration.hpp"
#include "common/logging/logger.hpp"
#include "common/maya/factory.hpp"
#include "common/maya/toonRenderer.hpp"
#include "common/util/definitions.hpp"
#include "common/util/math.hpp"
#include "common/util/maya.hpp"

namespace bookmark_squidget {
namespace apps {

namespace {
using AttrValueMap = common::util::AttrValueMap;
using PenPoints = common::util::PenPoints;
using MayaAttributeEditor = common::maya::MayaAttributeEditor;
} // namespace
/**
 * BookmarkSquidgetTool does the following:
 * - Deforming object from attribute
 * - From a given stroke, look through all attributes that would deform
 * - Selects object attributes from stroke points
 * - deforms attributes.
 *
 * Given a set of deformable objects Mi with attributes Ai and a stroke S_c, and
 * generating toon curve T_j from M_i*A_i, find an attribute config A'_i that
 * generates toon curves T'_j that minimizes curve_dist(S_c, t'_j). In other
 * words, find an attribute matrix that renders toon curves that's closest to
 * our stroke. Have a regularization of A'_i - A_i (possibly a L1
 * regularization).
 *
 * Things to address:
 * - Create toon curve renders T_j for MObject M_i
 * - Compare curve distance of T_j to S_c.
 *      - Compare in viewport domain.
 * - Find creation and transform attributes.
 *      - Question: how do we recognize creation vs transform attributes?  Maybe
 * there's a hierarchy.
 * - Optimize iteratively to find deformation.
 *
 *
 * squidget
 */

class BookmarkSquidgetTool : public MPxToolCommand {
public:
  // S P C W L o
  struct RenderMatrices {
    MMatrix localMatrix;
    MMatrix worldMatrix;
    MMatrix cameraMatrix;
    MMatrix projMatrix;
    MMatrix screenMatrix;
    MMatrix screenTransformMatrix;
  };

  struct CurveDistanceItem {
    double squared_norm;
    double matrix_deviation;
    RenderMatrices renderMats;
    MDagPath geomPath;
    MPointArray shadowCurve;
    double score;
    MObject curveObj;

    CurveDistanceItem() { score = -1; }

    CurveDistanceItem(double squared_norm, double matrix_deviation,
                      RenderMatrices renderMats, MDagPath geomPath,
                      MPointArray shadowCurve)
        : squared_norm(squared_norm), matrix_deviation(matrix_deviation),
          renderMats(renderMats), geomPath(geomPath), shadowCurve(shadowCurve) {
      double alpha = 0.65;
      score = (1 - alpha) * 1 / squared_norm + alpha * 1 / matrix_deviation;
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

  BookmarkSquidgetTool();
  ~BookmarkSquidgetTool() override;
  static void *creator();

  MStatus doIt(const MArgList &args) override;
  MStatus cancel() override;

  MStatus redoIt() override;
  MStatus undoIt() override;
  // bool isUndoable() const override;
  MStatus finalize() override;
  static MSyntax newSyntax();

  void ReadPoint(MPoint pt, MVector offsetVec, bool isCtrl);

  MStatus getControlPt(MPoint &controlPt);
  MStatus parseArgs(const MArgList &args);
  MStatus predoIt();

  void setPenControlPoints(MPointArray arr);
  void queryClosestObjectOnly();
  // void queryClosestObjectOnly(MPointArray strokePoints2D);
  // void queryClosestObject(MPointArray strokePoints2D);

  MPoint controlPoint;
  MDagPath controlPath;
  RenderMatrices controlMatrices;
  static CurveDistanceItem toDeformItem;

private:
  void deformClosestObject();

  MObjectArray getImplicitSquidgetCurves(MStringArray inputGeom,
                                         bool worldSpace);
  RenderMatrices getRenderMatrices(MObject geomObj);

  MStringArray getGeomsAroundStroke();
  MPointArray findShadowCurve(MPointArray strokePoints2D, MPointArray cvs,
                              RenderMatrices renderMats, double noise);

  // Eigen::Matrix4d leastSquares2(M4xX objectToWorld, M4xX screenToWorld);
  // Eigen::Matrix4d leastSquares(MPointArray shadowCurve, RenderMatrices
  // renderMats);
  Eigen::Matrix4d leastSquares(MPointArray initialPoints,
                               MPointArray targetPoints,
                               RenderMatrices renderMats);
  Eigen::Matrix3d leastSquaresTemp(M3xX x, M3xX b);

  double calculateShadowCurveError(RenderMatrices renderMats, MPointArray cvs);
  double matrixDeviation(MMatrix a, MMatrix b);
  bool crossedOver(MMatrix init, MMatrix post);
  bool overRotate(MMatrix init, MMatrix post);

  PenPoints penControlPoints; // screen Space
  // PenPoints queryControlPoints; // screen Space

  // MDagPath toDeformPath;
  // MPointArray toDeformShadowCurve;

  MTransformationMatrix oldMatrix;
  MVector oldDragDifference;

  // IMPORTANT: This is the main selection group used for the tool.
  int mode = 2; // study_0
                // 1; // study_1
                // 2; // study_2
};

} // namespace apps
} // namespace bookmark_squidget

#endif