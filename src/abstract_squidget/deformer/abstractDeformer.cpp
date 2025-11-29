#include "abstract_squidget/deformer/abstractDeformer.hpp"

#include <queue>

#include <igl/AABB.h>
#include <igl/cotmatrix.h>

#include <maya/M3dView.h>
#include <maya/MFnMesh.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnTransform.h>
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>

#include "abstract_squidget/deformer/deformerCurves.hpp"
#include "abstract_squidget/deformer/errorCalculation.hpp"
#include "abstract_squidget/deformer/method/deform.hpp"
#include "abstract_squidget/deformer/optimization.hpp"
#include "abstract_squidget/deformer/renderMatrices.hpp"
#include "common/util/maya.hpp"

#include "common/maya/factory.hpp"

/**
 * @brief NOTES
 * - ImplicitSquidgetCurves: We don't need the toon renderer to get these curves
 * but rather we can just use the densely sampled points.
 * - Shadow Curve = Curve to Curve matching:
 *    - Given a curve C, find C' that matches Q the best.  We can say that C has
 *      points [c1 ... cm], so there are O(m^2) possible subcurves.
 *    - 3 Energies:
 *      - Spatial Term: Best fit rigid or affine transformation
 *      - Alignment Term: Curvature + stretching/shrinking
 *      - shape energy: point to point distance
 *    - How to optimize for all energies.
 *
 */

namespace abstract_squidget {
namespace deformer {
QueryClosestObjectResult queryPreselectObject(MString object,
                                              MPointArray strokePoints2D,
                                              MPointArray shadowCurve,
                                              EditFlags flags) {
  MDagPath geomPath = common::util::getDagPathFromName(object);
  RenderMatrices renderMats = getRenderMatrices(geomPath.node());
  RenderMatrices initRenderMats = renderMats;

  renderMats =
      optimizeNewMatrices(shadowCurve, strokePoints2D, renderMats, flags);

  QueryClosestObjectResult result;
  result.object = geomPath.fullPathName();
  result.preShadowCurve = shadowCurve;

  MPointArray postShadowCurve;
  for (unsigned int i = 0; i < shadowCurve.length(); ++i) {
    MPoint p = renderMats.localMatrix * shadowCurve[i];
    p.cartesianize();
    postShadowCurve.append(p);
  }
  result.postShadowCurve = postShadowCurve;
  result.renderMats = renderMats;
  result.initRenderMats = initRenderMats;
  return result;
}

// =============================================================================
//                            Query and Rigid Transform
// =============================================================================
QueryClosestObjectResult
queryClosestObject(MStringArray candidateObjs, MPointArray curvePoints2D,
                   EditFlags flags,
                   std::shared_ptr<common::logging::Logger> logger) {
  /**
   * @brief Query closest object to the 2D stroke points.
   * Closest object is defined by the following:
   * - Distance from stroke points to object's shadow curve.
   * - Shape difference of object's shadow curve to stroke points.
   */
  MStringArray outputGeoms = candidateObjs; // Select closest Geometry

  std::priority_queue<CurveDistanceItem, std::vector<CurveDistanceItem>,
                      CompareCurveDistanceItem>
      pq;

  for (size_t i = 0; i < outputGeoms.length(); ++i) {
    MString outputGeom = outputGeoms[i];
    MDagPath geomPath = common::util::getDagPathFromName(outputGeom);
    MString inputGeoms[] = {outputGeom};
    MObjectArray abstractSquidgets = getImplicitSquidgetCurves(
        MStringArray(inputGeoms, 1), false, flags.isRigCurves);

    // For each squidget curve obj from the output geometry
    for (MObject squidgetObj : abstractSquidgets) {
      RenderMatrices renderMats = getRenderMatrices(geomPath.node());
      RenderMatrices initRenderMats = renderMats;
      MPointArray squidgetObjCVs = getObjSpaceCVs(squidgetObj, renderMats);
      if (squidgetObjCVs.length() < 2) {
        cout << "SHORT:" + geomPath.fullPathName() + "squidgetObjCVs.length(): "
             << squidgetObjCVs.length() << endl;
        continue;
      }

      MPointArray initShadowObjCVs =
          findShadowCurve(curvePoints2D, squidgetObjCVs, renderMats, 0);

      MPointArray shadowObjCVsScreenSpace;
      method::convertPointsToScreenSpace(initShadowObjCVs, renderMats,
                                         shadowObjCVsScreenSpace);

      renderMats = optimizeNewMatrices(initShadowObjCVs, curvePoints2D,
                                       renderMats, flags);

      // Project into screen space
      // I don't use worldToScreen because i don't remember why, but this is
      // more accurate.
      double squared_norm =
          calculateShadowCurveError(renderMats, squidgetObjCVs, curvePoints2D);
      double matrix_distance =
          matrixDeviation(initRenderMats.localMatrix, renderMats.localMatrix);

      CurveDistanceItem item(squared_norm, matrix_distance, renderMats,
                             initRenderMats, geomPath, squidgetObjCVs,
                             initShadowObjCVs);
      pq.push(item);

      MObject sCur = common::maya::CreateNurbsCurve(shadowObjCVsScreenSpace, 0);
      MDagPath curvePath = common::util::getNurbsCurveDagPathFromObject(sCur);
      logger->logCurveDistScore(curvePath.fullPathName().asChar(), item.score,
                                squared_norm, matrix_distance);
      common::util::deleteCurve(sCur);
      ;
    }

    // Delete all the squidgets from calculations
    for (MObject squidgetObj : abstractSquidgets) {
      common::util::deleteCurve(squidgetObj);
      ;
    }
  }

  QueryClosestObjectResult result;
  if (pq.size() == 0) { // cout << "pq is empty" << endl;
    logger->log("No objects to deform");
    result.object = MString("");
    result.preShadowCurve = MPointArray();
    result.postShadowCurve = MPointArray();
    result.renderMats = RenderMatrices();
  } else {

    CurveDistanceItem toDeformItem = pq.top();
    MGlobal::selectByName(toDeformItem.geomPath.fullPathName(),
                          MGlobal::kReplaceList);
    result.object = toDeformItem.geomPath.fullPathName();
    result.preShadowCurve = toDeformItem.shadowCurve;

    // // Convert preShadowCurve to world space
    // MPointArray preWorld;
    // MPointArray postWorld;
    // for (unsigned int i = 0; i < toDeformItem.shadowCurve.length(); ++i) {
    //   MPoint pr = toDeformItem.initRenderMats.worldMatrix *
    //   toDeformItem.initRenderMats.localMatrix * toDeformItem.shadowCurve[i];
    //   MPoint po = toDeformItem.renderMats.worldMatrix *
    //   toDeformItem.renderMats.localMatrix * toDeformItem.shadowCurve[i];
    //   preWorld.append(pr);
    //   postWorld.append(po);
    // }
    // common::maya::CreateNurbsCurve(preWorld, 0);
    // common::maya::CreateNurbsCurve(postWorld, 0);

    MPointArray postShadowCurve;
    for (unsigned int i = 0; i < toDeformItem.shadowCurve.length(); ++i) {
      MPoint p =
          toDeformItem.renderMats.localMatrix * toDeformItem.shadowCurve[i];
      p.cartesianize();
      postShadowCurve.append(p);
    }
    result.postShadowCurve = postShadowCurve;
    result.renderMats = toDeformItem.renderMats;
    result.initRenderMats = toDeformItem.initRenderMats;
  }
  return result;
}

DeformObjectResult transformObject(MString object, MPointArray curvePoints2D,
                                   RenderMatrices renderMats, EditFlags flags) {
  MDagPath geomPath = common::util::getDagPathFromName(object);

  // Apply the new matrix
  MFnTransform fnMinTransform;
  fnMinTransform.setObject(geomPath);
  MTransformationMatrix oldTransform = fnMinTransform.transformation();
  MTransformationMatrix newTransform(renderMats.localMatrix.transpose());
  if (flags.translateOnly) {
    MQuaternion oldQ = oldTransform.rotation();
    newTransform.setRotationQuaternion(oldQ.x, oldQ.y, oldQ.z, oldQ.w);
  }
  double scale[3] = {1, 1, 1};
  double sheer[3] = {0, 0, 0}; // Maintain Sheer
  oldTransform.getScale(scale, MSpace::kPostTransform);
  oldTransform.getShear(sheer, MSpace::kPostTransform);
  newTransform.setScale(scale, MSpace::kPostTransform);
  newTransform.setShear(sheer, MSpace::kPostTransform);

  fnMinTransform.set(newTransform.asMatrix());

  MStatus s;
  MPoint controlPoint =
      fnMinTransform.getTranslation(MSpace::kPostTransform, &s);
  MPoint controlScreenPt = (renderMats.screenMatrix * renderMats.projMatrix *
                            renderMats.cameraMatrix * renderMats.worldMatrix) *
                           controlPoint;
  controlScreenPt.cartesianize();

  MVector oldDragDifference =
      curvePoints2D[curvePoints2D.length() - 1] -
      controlScreenPt; // Save old drag difference for screen space holding.
  oldDragDifference.z = 0;

  DeformObjectResult result;
  result.deformType = RIGID_TRANSFORM;
  result.oldTransform = oldTransform;
  result.offset = oldDragDifference;
  return result;
}

// =============================================================================
//                            Vertex Deformation
// =============================================================================
DeformObjectResult deformObject(MString object, MPointArray _curvePoints2D,
                                MPointArray shadowCurve,
                                RenderMatrices renderMats, EditFlags flags) {
  /*
    Convert shadowCurve to screen space.
    Get object vertices closest to shadowCurve2D.
      - Get connected interval of points?
    Morph  vertices in 2D to match the curvePoints2D.
    Apply transformation to those vertices in object space
  */

  MPointArray curvePoints2D;
  // smooth curvePoints2D once
  for (unsigned int i = 1; i < _curvePoints2D.length() - 1; ++i) {
    MPoint p =
        (_curvePoints2D[i - 1] + _curvePoints2D[i] + _curvePoints2D[i + 1]) /
        3.0;
    curvePoints2D.append(p);
  }

  // 0) Get the mesh vertices in local space
  MDagPath geomPath = common::util::getDagPathFromName(object);
  MPointArray meshVertices;
  MFnMesh fnMesh(geomPath);
  fnMesh.getPoints(meshVertices, MSpace::kObject);
  int n = meshVertices.length();

  // 1) Convert mesh vertices to screen space
  MPointArray meshVerticesScreenSpace;
  method::convertPointsToScreenSpace(meshVertices, renderMats,
                                     meshVerticesScreenSpace);
  MPointArray shadowCurveScreenSpace;
  method::convertPointsToScreenSpace(shadowCurve, renderMats,
                                     shadowCurveScreenSpace);

  // 2) Find the closest vertices to the shadow curve by the interval
  int start_interval;
  int end_interval;
  int interval_direction;
  method::findInterval(meshVerticesScreenSpace, shadowCurveScreenSpace,
                       start_interval, end_interval, interval_direction);

  // 2.5) Get the interval curve
  MPointArray meshCurveIntervalScreen;
  int curr_i = start_interval;
  do {
    meshCurveIntervalScreen.append(meshVerticesScreenSpace[curr_i]);
    curr_i = (curr_i + interval_direction + n) % n;
  } while (curr_i != end_interval);

  // 3) Parameterize mesh vertices from [0, 1]
  // X and Y have same orientation
  Eigen::MatrixX2d X =
      method::convertMPointArrayToMatrixX2d(meshCurveIntervalScreen);
  Eigen::MatrixX2d Y = method::convertMPointArrayToMatrixX2d(curvePoints2D);
  Eigen::VectorXd X_t = method::parameterize_curve(X);
  Eigen::VectorXd Y_t = method::parameterize_curve(Y);

  // 3.5) Set the falloff distance so there's not overlap
  // softSelect -e -softSelectDistance 2.347418;
  std::ostringstream ss;
  ss << "softSelect -q -softSelectDistance";
  MString res = MGlobal::executeCommandStringResult(ss.str().c_str());
  double old_falloff_distance = res.asDouble();

  // 3.6) Set the falloff distance to half the distance between the first and
  // last point
  MPoint first_point, last_point;
  fnMesh.getPoint(start_interval, first_point, MSpace::kObject);
  fnMesh.getPoint(end_interval, last_point, MSpace::kObject);
  double temp_falloff_distance = (last_point - first_point).length() / 2.0;
  ss.str("");
  ss << "softSelect -e -softSelectDistance " << temp_falloff_distance << ";";
  MGlobal::executeCommand(ss.str().c_str());

  // 4) Transform all vertices of X to fit Y
  // cout << "handle positions" << endl;
  Eigen::MatrixX3d X_after(X.rows(), 3);
  MPointArray handle_positions;
  for (int i = 0; i < X.rows(); ++i) {
    curr_i = (start_interval + i * interval_direction + n) % n;
    // cout << "   " << curr_i << " ";

    MPoint p;
    fnMesh.getPoint(curr_i, p, MSpace::kObject);
    MPoint p2d;
    method::convertPointToScreenSpace(p, renderMats, p2d);

    // Find the point along curve using param
    double param = X_t(i);
    Eigen::Vector2d X0(p2d.x, p2d.y);
    Eigen::Vector2d X1 = method::findPointOnCurve(Y, Y_t, param);
    // cout << X0 << " " << X1 << endl;

    // x0 = S P C W L m
    // D S P C W L = S P C W L B m
    // D (S P C W L) = (S P C W L) B
    // D A = A B
    // A^-1 D A = B
    MVector d(X1.x() - X0.x(), X1.y() - X0.y(), 0);
    MMatrix dAsMatrix = MMatrix::identity;
    dAsMatrix(0, 3) = d.x;
    dAsMatrix(1, 3) = d.y;
    dAsMatrix(2, 3) = d.z;
    MMatrix A = renderMats.screenMatrix * renderMats.projMatrix *
                renderMats.cameraMatrix * renderMats.worldMatrix *
                renderMats.localMatrix;
    MMatrix B = A.inverse() * dAsMatrix * A;

    // Get handle positions
    MPoint p2 = B * p;
    handle_positions.append(p2);
    X_after.row(i) = Eigen::Vector3d(p2.x, p2.y, p2.z);

    // Maybe set the falloff distance
    if (i == 0 || i == X.rows() - 1) {
      MPoint diff = renderMats.localMatrix * (p2 - p);

      ss.str("");
      ss << "select -r " << object << ".vtx[" << curr_i << "];";
      ss << "move -r " << diff.x << " " << diff.y << " " << diff.z << ";";
      ss << endl;
      // cout << ss.str() << endl;
      MGlobal::executeCommand(ss.str().c_str());
    }
  }

  MPointArray postShadowCurve;
  for (int i = 0; i < Y_t.rows(); ++i) {
    double param = Y_t(i);
    Eigen::Vector3d x = method::findPointOnCurve(X_after, X_t, param);
    postShadowCurve.append(MPoint(x.x(), x.y(), x.z()));
  }

  // deformFalloff(fnMesh, handle_positions, 0, start_interval,
  // interval_direction,
  //               renderMats);
  // deformFalloff(fnMesh, handle_positions, X.rows() - 1, start_interval,
  //               interval_direction, renderMats);

  // Move interior points to prevent local deformation.
  for (int i = 0; i < X.rows(); ++i) {
    curr_i = (start_interval + i * interval_direction + meshVertices.length()) %
             meshVertices.length();
    MPoint p2 = handle_positions[i];
    fnMesh.setPoint(curr_i, p2, MSpace::kObject);
  }

  // ss.str("");
  // ss << "softSelect -e -softSelectDistance " << old_falloff_distance << ";";
  // MGlobal::executeCommand(ss.str().c_str());

  DeformObjectResult result;
  result.deformType = CURVE_DEFORM;
  result.oldVertices = meshVertices;
  result.postShadowCurve = postShadowCurve;
  return result;
}

void deformFalloff(MFnMesh &fnMesh, MPointArray handle_positions,
                   int start_stroke_i, int start_curve_i,
                   int interval_direction, RenderMatrices renderMats) {
  int n = fnMesh.numVertices();

  int a_i, b_i;
  int falloff_direction;
  if (start_stroke_i == 0) {
    a_i = 0;
    b_i = 1;
    falloff_direction = (interval_direction > 0) ? -1 : 1;
  } else {
    a_i = start_stroke_i;
    b_i = start_stroke_i - 1;
    falloff_direction = interval_direction;
  }

  // Get the first edge of the curve in screen space
  MPoint XA, XB;
  fnMesh.getPoint((start_curve_i + a_i * interval_direction + n) % n, XA,
                  MSpace::kObject);
  fnMesh.getPoint((start_curve_i + b_i * interval_direction + n) % n, XB,
                  MSpace::kObject);
  MPoint XA2D, XB2D;
  method::convertPointToScreenSpace(XA, renderMats, XA2D);
  method::convertPointToScreenSpace(XB, renderMats, XB2D);

  // Get the first edge of the stroke in screen space
  MPoint YA = handle_positions[a_i], YB = handle_positions[b_i];
  MPoint YA2D, YB2D;
  method::convertPointToScreenSpace(YA, renderMats, YA2D);
  method::convertPointToScreenSpace(YB, renderMats, YB2D);

  // Rotation in screen space
  MVector X_edge_2d = XB2D - XA2D;
  MVector Y_edge_2d = YB2D - YA2D;

  double cos_theta = X_edge_2d.normal() * Y_edge_2d.normal();
  double theta = acos(cos_theta);

  double cross_product = X_edge_2d.x * Y_edge_2d.y - X_edge_2d.y * Y_edge_2d.x;
  if (cross_product < 0) {
    theta = -theta;
  }

  MMatrix R = MMatrix::identity;
  R(0, 0) = cos(theta);
  R(0, 1) = -sin(theta);
  R(1, 0) = sin(theta);
  R(1, 1) = cos(theta);

  // Translation in screen space
  MVector translation = YA2D - R * XA2D;

  // Apply fall off to vertices beyond edge
  int falloff_points = 200;
  for (int i = 0; i < falloff_points; ++i) {
    int curr_i = (start_curve_i + (a_i + i) * falloff_direction + n) % n;

    // double alpha = 1.0 - sigmoid((double)i / falloff_points);
    double alpha = 1.0 - quintic((double)i / falloff_points);
    MMatrix dAsMatrix =
        screenTransformationMatrix(alpha * translation, alpha * theta);
    MMatrix A = renderMats.screenMatrix * renderMats.projMatrix *
                renderMats.cameraMatrix * renderMats.worldMatrix *
                renderMats.localMatrix;
    MMatrix B = A.inverse() * dAsMatrix * A;

    MPoint p;
    fnMesh.getPoint(curr_i, p, MSpace::kObject);
    MPoint p2 = B * p;
    fnMesh.setPoint(curr_i, p2, MSpace::kObject);
  }
}

MMatrix screenTransformationMatrix(MVector t, double angle) {
  MMatrix dAsMatrix = MMatrix::identity;
  dAsMatrix(0, 3) = t.x;
  dAsMatrix(1, 3) = t.y;
  dAsMatrix(2, 3) = t.z;
  dAsMatrix(0, 0) = cos(angle);
  dAsMatrix(0, 1) = -sin(angle);
  dAsMatrix(1, 0) = sin(angle);
  dAsMatrix(1, 1) = cos(angle);
  return dAsMatrix;
}

double sigmoid(double x) { return 1.0 / (1.0 + exp(-10 * (x - 0.5))); }
double quintic(double x) {
  return 6 * pow(x, 5) - 15 * pow(x, 4) + 10 * pow(x, 3);
}

// =============================================================================
//                            Other
// =============================================================================
MVector livePoint(MString liveObject, MPoint mousePt, MVector oldDragDifference,
                  bool isCtrl) {
  // Assume controlObj exists
  // cout << "abstractDeformer::livePoint()" << endl;

  MDagPath geomPath = common::util::getDagPathFromName(liveObject);
  MFnTransform fnTransform(geomPath);
  RenderMatrices controlMatrices = getRenderMatrices(geomPath.node());

  MStatus s;
  MPoint controlPoint = fnTransform.getTranslation(MSpace::kWorld, &s);
  MPoint controlScreenPt =
      (controlMatrices.screenMatrix * controlMatrices.projMatrix *
       controlMatrices.cameraMatrix) *
      controlPoint;

  M4x4 screenATrans = M4x4::Identity();
  M4x4 screenARot = M4x4::Identity();

  M4x4 T = M4x4::Identity();
  M4x4 R = M4x4::Identity();
  if (isCtrl) {
    MVector difference = mousePt - controlScreenPt;
    difference.z = 0;
    // Check if difference is rotated clockwise to oldDragDifference
    double theta = difference.angle(oldDragDifference);
    bool isClockwise = difference.x * oldDragDifference.y -
                           difference.y * oldDragDifference.x >
                       0;
    if (isClockwise) {
      theta = -theta;
    }

    oldDragDifference = difference;

    screenARot(0, 0) = cos(theta);
    screenARot(0, 1) = -sin(theta);
    screenARot(1, 0) = sin(theta);
    screenARot(1, 1) = cos(theta);

    // A S' S P C W = S' S P C (T W R)
    // A S' S P C W = S' S P C (T W R)
    // A X = X R
    // X^-1 A X = R
    M4x4 A = screenARot;
    M4x4 screenTranslate = M4x4::Identity();
    screenTranslate(0, 3) = -controlMatrices.screenMatrix(0, 3);
    screenTranslate(1, 3) = -controlMatrices.screenMatrix(1, 3);

    M4x4 localTranslate = M4x4::Identity();
    localTranslate(0, 3) = -controlMatrices.localMatrix(0, 3);
    localTranslate(1, 3) = -controlMatrices.localMatrix(1, 3);
    localTranslate(2, 3) = -controlMatrices.localMatrix(2, 3);

    M4x4 X = screenTranslate *
             common::util::mMatrixToMatrix(controlMatrices.screenMatrix) *
             common::util::mMatrixToMatrix(controlMatrices.projMatrix) *
             common::util::mMatrixToMatrix(controlMatrices.cameraMatrix) *
             common::util::mMatrixToMatrix(controlMatrices.worldMatrix) *
             localTranslate;
    R = X.inverse() * A * X;

    // A S' S P C = S' S P C -T R T
    // X =  common::util::mMatrixToMatrix(controlMatrices.screenMatrix) *
    //       common::util::mMatrixToMatrix(controlMatrices.projMatrix) *
    //       common::util::mMatrixToMatrix(controlMatrices.cameraMatrix);
    // A S' X = S' X -T R T
    // R = X.inverse() * A * X;
    MVector oldTrans;
    oldTrans.x = controlMatrices.localMatrix(0, 3);
    oldTrans.y = controlMatrices.localMatrix(1, 3);
    oldTrans.z = controlMatrices.localMatrix(2, 3);

    controlMatrices.localMatrix =
        common::util::matrixToMMatrix(R) * controlMatrices.localMatrix;
    controlMatrices.localMatrix(0, 3) = oldTrans.x;
    controlMatrices.localMatrix(1, 3) = oldTrans.y;
    controlMatrices.localMatrix(2, 3) = oldTrans.z;

  } else {
    MVector newDragdifference = mousePt - controlScreenPt;
    MVector toTranslate = newDragdifference - oldDragDifference;

    // Translate only by mouseObjectOffsetPoint and currentObjectPoint
    // MPoint offsetMousePt = mousePt + offsetVec;
    screenATrans(0, 3) = toTranslate.x;
    screenATrans(1, 3) = toTranslate.y;

    // A S P C W L = S P C W (T L R)
    // A S P C W L = S P C W T L
    // A X = X T
    // X^-1 A X = T
    M4x4 A = screenATrans;
    M4x4 X = common::util::mMatrixToMatrix(controlMatrices.screenMatrix) *
             common::util::mMatrixToMatrix(controlMatrices.projMatrix) *
             common::util::mMatrixToMatrix(controlMatrices.cameraMatrix) *
             common::util::mMatrixToMatrix(controlMatrices.worldMatrix);
    T = X.inverse() * A * X;
    controlMatrices.localMatrix =
        common::util::matrixToMMatrix(T) * controlMatrices.localMatrix;
  }

  fnTransform.set(MTransformationMatrix(controlMatrices.localMatrix.transpose())
                      .asMatrix());
  double scale[3] = {1, 1, 1};
  fnTransform.setScale(scale);
  double sheer[3] = {0, 0, 0};
  fnTransform.setShear(sheer);
  return oldDragDifference;
}

MPointArray getObjSpaceCVs(MObject nurbsCurveObj, RenderMatrices renderMats) {
  /**
   * @brief Return the control vertices in object space of a NURBS curve
   * object.
   */
  // cout << "abstractDeformer::getObjSpaceCVs()" << endl;
  MDagPath objPath =
      common::util::getNurbsCurveDagPathFromObject(nurbsCurveObj);
  MFnNurbsCurve curveFn(objPath);

  MPointArray worldSpaceCVs;
  curveFn.getCVs(worldSpaceCVs);

  // Converts CVs to object space to match object matrix
  MPointArray objSpaceCVs;
  for (MPoint p : worldSpaceCVs) {
    MPoint objp =
        (renderMats.worldMatrix * renderMats.localMatrix).inverse() * p;
    objSpaceCVs.append(objp);
  }
  return objSpaceCVs;
}
} // namespace deformer
} // namespace abstract_squidget
