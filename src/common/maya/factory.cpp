#include "common/maya/factory.hpp"

namespace common {
namespace maya {
MObject CreateNurbsCurve(MPointArray stroke, int deg, MObject &parent) {
  MStatus status;
  int ncvs = stroke.length();
  int order = deg + 1;
  int nknots = ncvs + order;

  MPointArray curvePoints;
  MDoubleArray knotSequence;

  for (unsigned j = 0; j < stroke.length(); ++j) {
    curvePoints.append(MPoint(stroke[j]));
  }

  // Create knot sequence
  for (int j = 0; j < order - 1; ++j) {
    knotSequence.append(0);
  }
  for (int j = 0; j < nknots - 2 * (order - 1); ++j) {
    // double t = (double) j / (ncvs - order);
    double t = (double)j;
    knotSequence.append(t);
  }
  // for (int j = 0; j < order-1; ++j) { knotSequence.append(1); }
  for (int j = 0; j < order - 1; ++j) {
    knotSequence.append(ncvs - order);
  }

  bool create2D = false;
  bool createRational = false;
  bool uniformParam = true;
  MFnNurbsCurve curveFn;
  MObject curve = curveFn.createWithEditPoints(
      curvePoints, order, MFnNurbsCurve::kOpen, create2D, createRational,
      uniformParam, MObject::kNullObj, &status);

  // Set new parent object
  MFnDagNode fnDagNode;
  fnDagNode.setObject(parent);
  fnDagNode.addChild(curve);

  return curve;
}

} // namespace maya
} // namespace common
