#include "abstract_squidget/deformer/errorCalculation.hpp"

#include <maya/MFnNurbsCurve.h>

#include <abstract_squidget/deformer/method/deform.hpp>
#include <common/maya/toonRenderer.hpp>
#include <common/util/maya.hpp>

// =============================================================================
//                           Error Calculation
// =============================================================================
namespace abstract_squidget {
namespace deformer {
double calculateShadowCurveError(RenderMatrices renderMats, MPointArray objCvs,
                                 MPointArray curvePoints2D) {
  // Convert cvs to screen space
  MPointArray objScreenCVs;
  method::convertPointsToScreenSpace(objCvs, renderMats, objScreenCVs);

  // Find closest point of screenCVs on curve to strokePoints
  // Find params of closest points.  Map params to curve.
  MObject objCurve = common::maya::CreateNurbsCurve(objScreenCVs, 0);
  MDagPath objCurvePath =
      common::util::getNurbsCurveDagPathFromObject(objCurve);
  MFnNurbsCurve objCurveFn(objCurvePath);

  // curveFn.setCVs(screenCVs);
  double totalDistance = 0;
  MPointArray shadowCVs;
  for (MPoint p : curvePoints2D) {
    MPoint shadowPoint = objCurveFn.closestPoint(p);
    shadowCVs.append(shadowPoint);
    totalDistance += (shadowPoint - p).length();
  }



  common::util::deleteCurve(objCurve);
  return totalDistance / curvePoints2D.length() / 2;
}

double matrixDeviation(MMatrix a, MMatrix b) {
  MVector aT(a[0][3], a[1][3], a[2][3]);
  MVector bT(b[0][3], b[1][3], b[2][3]);

  // return pow((aT - bT).length(), 2);
  return (aT - bT).length();
}

bool crossedOver(MPointArray curvePoints2D, MMatrix a, MMatrix b) {
  bool crossesOver = false;
  MPoint aT(a[0][3], a[1][3], a[2][3]);
  MPoint bT(b[0][3], b[1][3], b[2][3]);

  M3dView view = M3dView::active3dView();
  short ax, ay, bx, by;
  view.worldToView(aT, ax, ay);
  view.worldToView(bT, bx, by);
  MPoint aScreen(ax, ay, 0);
  MPoint bScreen(bx, by, 0);

  MPoint strokeMidpoint(0, 0, 0);
  for (MPoint p : curvePoints2D) {
    strokeMidpoint += p;
  }
  strokeMidpoint = strokeMidpoint / curvePoints2D.length();

  MVector AP = aScreen - strokeMidpoint;
  MVector BP = bScreen - strokeMidpoint;

  double angle = AP.angle(BP);
  // cout << "angle: " << angle << " " << aScreen << bScreen << strokeMidpoint
  // << endl;
  return angle > M_PI / 2.;
}

bool overRotate(MMatrix a, MMatrix b) {
  // checks if rotation of a to b is greater than 90 degrees
  MTransformationMatrix aTrans(a);
  MTransformationMatrix bTrans(b);

  MVector aRot = aTrans.eulerRotation().asVector();
  MVector bRot = bTrans.eulerRotation().asVector();

  double angle = aRot.angle(bRot);
  // //cout << "angle: " << angle << endl;
  return angle > M_PI / 2.;
}

} // namespace deformer
} // namespace abstract_squidget