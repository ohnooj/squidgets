#include "bookmark_squidget/canvas/interface/elementContainer.hpp"

// Assume all points coming in are in local space.
namespace bookmark_squidget {
namespace canvas {
namespace interface {
ElementContainer::ElementContainer() {
  canvasShape.extendToShape();
}

ElementContainer::ElementContainer(MDagPath canvasPath) {
  canvasShape = canvasPath;
  canvasShape.extendToShape();
}

ElementContainer::~ElementContainer() {
  MObject o = pathShape.transform();
  MGlobal::deleteNode(o);

  o = ghostShape.transform();
  MGlobal::deleteNode(o);
}

// =============================================================================
// Add Step Elements
// =============================================================================
void ElementContainer::AddStepElement(
    std::shared_ptr<StepElement> stepElement) {
  m_stepElements.push_back(stepElement);
  updateInterface();
  for (auto stepElement : m_stepElements) {
    MFnDagNode fnDagNode(stepElement->curve);
  }
}

// =============================================================================
// Delete Step Elements
// =============================================================================
// void ElementContainer::DeleteStepElement(int index) {
//   // m_stepElements.
// cout << "Deleting: " << index << endl;
//   m_stepElements.erase(m_stepElements.begin() + index);
//   updateInterface();
// }

double ElementContainer::UpdateValueAlongPath(MPoint localPt) {
  if (!pathShape.isValid()) {
    return -1;
  }

  double t = common::util::getDistanceAlongCurve(localPt, pathShape.node());
  t = t / (m_stepElements.size() - 1);

  MFnDagNode fnDagNode(pathShape);
  MPlug tPlug = fnDagNode.findPlug(squidget_t, true);
  double du;
  tPlug.getValue(du);
  tPlug.setDouble(t);

  return t;
}

// =============================================================================
// Update Interface
// =============================================================================
void ElementContainer::updateInterface() {
  bool pathInstanceExists = pathShape.isValid();
  bool pathLengthValid = pathLongEnough();

  if (!pathInstanceExists && pathLongEnough()) {
    createInterface();
  } else if (pathInstanceExists && !pathLongEnough()) {
    deleteInterface();
  } else if (pathInstanceExists && pathLongEnough()) {
    deleteInterface();
    createInterface();
  }
}

// =============================================================================
// Create Interface
// =============================================================================
void ElementContainer::createInterface() {
  createPath();
  createTAttribute();
  createGhost();
  createElementBlendShape();
  connectSquidgetToBlendshape();
}

void ElementContainer::createPath() {
  // TODO: Possibly update the CV's instead of create/delete

  // Find points of Path from midpoints of Step curves
  MPointArray pathPoints;
  for (size_t i = 0; i < m_stepElements.size(); i++) {
    auto stepElement = m_stepElements[i];

    MDagPath path;
    MFnDagNode fnDagNode(stepElement->curve);
    fnDagNode.getPath(path);
    path.extendToShapeDirectlyBelow(0);

    MPoint mid = common::util::getPointAlongCurve(0.5, path.node());
    pathPoints.append(mid);
  }

  // Create path curve parented to canvasObj
  MObject canvasObj = canvasShape.transform();
  MObject curve = common::maya::CreateNurbsCurve(pathPoints, 0, canvasObj);

  // Save path transform
  MFnDagNode fnDagNode(curve);
  fnDagNode.getPath(pathShape);
  pathShape.extendToShapeDirectlyBelow(0);

  // Renaming Path
  MDGModifier modifier;
  modifier.renameNode(curve, "Path");
  modifier.doIt();
}

void ElementContainer::createTAttribute() {
  MFnNumericAttribute attr;
  squidget_t = attr.create("squidget_t", "squidget_t", MFnNumericData::kFloat);
  attr.setMin(0.0);
  attr.setMax(1.0);
  attr.setKeyable(true);
  attr.setChannelBox(true);
  attr.setWritable(true);
  attr.setReadable(true);

  MFnDagNode fnDagNode(pathShape);
  fnDagNode.addAttribute(squidget_t);
}

void ElementContainer::createGhost() {
  if (m_stepElements.size() < 2) {
    return;
  }

  auto baseStepElement = m_stepElements[0];
  MFnDagNode fnDagNode(baseStepElement->curve);

  MObject newNode = fnDagNode.duplicate();
  fnDagNode.setObject(newNode);
  std::ostringstream ss;
  ss << "rename " << fnDagNode.fullPathName() << " ";
  ss << "Ghost";
  MGlobal::executeCommand(ss.str().c_str());

  fnDagNode.getPath(ghostShape);
  ghostShape.extendToShape();
}

void ElementContainer::createElementBlendShape() {
  MDagPath dagPath;
  MFnDagNode fnDagNode;
  MObject baseShapeObj = ghostShape.node();

  MObject bsObj = deform.create(baseShapeObj);
  m_stepElementKey.clear();

  for (size_t i = 0; i < m_stepElements.size(); ++i) {
    auto se = m_stepElements[i];

    double weight = ((double)i) / (m_stepElements.size() - 1);
    weight = std::ceil(weight * 100.0) / 100.0;
    fnDagNode.setObject(se->curve);
    fnDagNode.getPath(dagPath);
    dagPath.extendToShape();

    deform.addTarget(baseShapeObj, 0, dagPath.node(), weight);
    m_stepElementKey.push_back(weight);
  }
}

void ElementContainer::connectSquidgetToBlendshape() {
  MFnDagNode fnDagNode(pathShape);
  MPlug tPlug = fnDagNode.findPlug(squidget_t, true);
  MPlug bsPlug = deform.findPlug(MString("weight"), true);
  MPlug wPlug = bsPlug.elementByPhysicalIndex(0);

  MDGModifier mod;
  MStatus stat = mod.connect(tPlug, wPlug);
  mod.doIt();
}

// =============================================================================
// Delete Interface
// =============================================================================
void ElementContainer::deleteInterface() {
  // cout << "deletePath()" << endl;
  MObject n = pathShape.transform();
  MGlobal::deleteNode(n);

  n = ghostShape.transform();
  MGlobal::deleteNode(n);
}

// =============================================================================
// Extra
// =============================================================================
double ElementContainer::queryStrokeValue(const MPointArray localPts,
                                          double &dist) {
  MFnDagNode fnDagNode;
  MDagPath path;

  fnDagNode.setObject(pathShape);
  MPlug tPlug = fnDagNode.findPlug(squidget_t, true);
  double oldWeight = tPlug.asDouble();

  // Get stroke as nurbs curve
  MObject canvasObj = canvasShape.transform();
  MObject curve = common::maya::CreateNurbsCurve(localPts, 0, canvasObj);
  curve = common::util::resampleCurve(curve, 30);
  fnDagNode.setObject(curve);
  fnDagNode.getPath(path);
  path.extendToShape();
  curve = path.node();

  // Upsidedown curve?  Who cares
  // Get Base object
  MObjectArray baseArr;
  deform.getBaseObjects(baseArr);
  MObject base = baseArr[0];

  // Getting blendShape input geometry and weights
  // blendShape1.inputTarget[0].inputTargetGroup[0].inputTargetItem[5000].inputGeomTarget
  MPlug inputTarget = deform.findPlug("inputTarget", false)[0];
  MPlug inputTargetGroup = inputTarget.child(0)[0];
  MPlug inputTargetItems = inputTargetGroup.child(0);

  // cout << "close curve time " << endl;
  double t = 0;
  double min_dist = std::numeric_limits<double>::max();
  MIntArray targetItemIndices;
  deform.targetItemIndexList(0, base, targetItemIndices);
  for (unsigned int i = 1; i < targetItemIndices.length(); ++i) {
    // Get two StepElements at a time.
    int b0_index = targetItemIndices[i - 1];
    double b0_t = (b0_index - 5000) / 1000.;
    MPlug b0_plug =
        inputTargetItems.elementByLogicalIndex(b0_index).child(0).source();
    MObject b0obj = b0_plug.node();
    fnDagNode.setObject(b0obj);

    Eigen::VectorXd b0 = common::util::nurbsToVector(b0obj);
    Eigen::VectorXd s0 = common::util::nurbsToVector(curve);
    b0 = translateToOrigin(b0);
    s0 = translateToOrigin(s0);

    Eigen::VectorXd s_hat = s0 - b0;

    int b1_index = targetItemIndices[i];
    double b1_t = (b1_index - 5000) / 1000.;
    MPlug b1_plug =
        inputTargetItems.elementByLogicalIndex(b1_index).child(0).source();
    MObject b1obj = b1_plug.node();
    fnDagNode.setObject(b1obj);

    Eigen::VectorXd b1 = common::util::nurbsToVector(b1obj);
    b1 = translateToOrigin(b1);
    Eigen::VectorXd db = b1 - b0;

    // w = (B^T * B)' * (B^T * s_hat)
    double denom = db.dot(db);
    double curr_t = 1. / denom * db.dot(s_hat);
    if (denom == 0) { // b1 and b0 are the same probably (single squidget)
      curr_t = 0;
    } else {
      curr_t = (b1_t - b0_t) * curr_t + b0_t;
    }

    tPlug.setDouble(curr_t);

    // Calculate distance
    Eigen::VectorXd f = common::util::nurbsToVector(ghostShape.node());
    double curr_dist = curveDistance(f, s0);

    if (curr_dist < min_dist) {
      min_dist = curr_dist;
      t = curr_t;
    }
  }

  t = t / (m_stepElements.size() - 1);
  t = std::max<double>(std::min<double>(t, 1), 0);

  tPlug.setDouble(oldWeight);
  MObject curveObj = path.transform();
  MGlobal::deleteNode(curveObj);
  dist = min_dist;
  return t;
}

Eigen::VectorXd ElementContainer::translateToOrigin(Eigen::VectorXd f) {
  Eigen::MatrixXd mf(3, f.size() / 3);
  for (int i = 0; i < mf.cols(); i++) {
    mf(0, i) = f(3 * i + 0);
    mf(1, i) = f(3 * i + 1);
    mf(2, i) = f(3 * i + 2);
  }

  // Calculate Average piont of f and s0
  auto mf_avg = mf.rowwise().mean();
  Eigen::MatrixXd origin_mf = mf.colwise() - mf_avg;

  Eigen::VectorXd origin_f(f.size());
  for (int i = 0; i < origin_mf.cols(); i++) {
    origin_f(3 * i + 0) = origin_mf(0, i);
    origin_f(3 * i + 1) = origin_mf(1, i);
    origin_f(3 * i + 2) = origin_mf(2, i);
  }

  return origin_f;
}

double ElementContainer::curveDistance(Eigen::VectorXd a, Eigen::VectorXd b) {
  Eigen::VectorXd origin_a = translateToOrigin(a);
  Eigen::VectorXd origin_b = translateToOrigin(b);
  Eigen::VectorXd diffs = origin_a - origin_b;
  double curr_dist = diffs.array().square().sum();

  // // reverse point order for one curve
  // Eigen::VectorXd rev_origin_a(origin_a.size());
  // for (int i = 0; i < rev_origin_a.size()/3; ++i) {
  //   int o = (rev_origin_a.size()/3)-1-i;
  //   rev_origin_a.block(3*i, 0, 3, 1) = origin_a.block(3*o, 0, 3, 1);
  // }
  // cout << "origin_a: " << endl << origin_a << endl;
  // cout << "rev_origin_a: " << endl << rev_origin_a << endl;

  // auto rev_diffs = origin_b - rev_origin_a;
  // double rev_curr_dist = diffs.array().square().sum();

  // curr_dist is min of curr_dist and rev_curr_dist
  // return std::min(curr_dist, rev_curr_dist);
  return curr_dist;
}

bool ElementContainer::setTValue(double t) {
  MFnDagNode fnDagNode(pathShape);
  MPlug tPlug = fnDagNode.findPlug(squidget_t, true);
  tPlug.setDouble(t);
  return true;
}

bool ElementContainer::pathLongEnough() { return m_stepElements.size() >= 2; }

int ElementContainer::CheckStrokeCrossing(MPointArray localStroke) {
  localStroke = common::util::swapPointArrayCoordinates(localStroke, 1, 2);
  for (size_t i = 0; i < m_stepElements.size(); i++) {
    auto stepElement = m_stepElements[i];
    MPointArray stepShape =
        common::util::swapPointArrayCoordinates(stepElement->shape, 1, 2);

    if (common::util::polylinesIntersect(localStroke, stepShape)) {
      return i;
    }
  }
  return -1;
}

double ElementContainer::CheckPathCrossing(MPointArray localStroke) {
  MPointArray stroke =
      common::util::swapPointArrayCoordinates(localStroke, 1, 2);

  MFnNurbsCurve fnCurve(pathShape);
  MPointArray pathPoints;
  fnCurve.getCVs(pathPoints);
  pathPoints = common::util::swapPointArrayCoordinates(pathPoints, 1, 2);

  MPoint closestPoint = localStroke[0];
  double minDist = std::numeric_limits<double>::max();
  double minParam = 1;

  for (size_t i = 0; i < localStroke.length(); i++) {
    double param;
    MPoint localPoint = localStroke[i];
    MPoint closestPoint = fnCurve.closestPoint(localPoint, &param);

    double dist = localPoint.distanceTo(closestPoint);

    if (dist < minDist) {
      minDist = dist;
      closestPoint = localPoint;
      minParam = param;
    }
  }

  if (minDist > 0.1) {
    return -1;
  }

  return minParam / (m_stepElements.size() - 1);
}
// double ElementContainer::CheckPathCrossing(MPointArray localStroke) {
//   // Return -1 if doesn't cross path, otherwise get crossing value along
//   curve. MPointArray stroke =
//   common::util::swapPointArrayCoordinates(localStroke, 1, 2);

//   MFnNurbsCurve fnCurve(pathShape);
//   MPointArray pathPoints;
//   fnCurve.getCVs(pathPoints);
//   pathPoints = common::util::swapPointArrayCoordinates(pathPoints, 1,
//   2); cout << "PathCrossing: " << pathPoints << endl;

//   for (size_t i = 0; i < m_stepElements.size(); i++) {
//     auto stepElement = m_stepElements[i];
//     MPointArray stepShape =
//     common::util::swapPointArrayCoordinates(stepElement->shape, 1,
//     2);

//     if (common::util::polylinesIntersect(stroke, stepShape)) {
//       return i;
//     }
//   }
//   for (size_t i = 1; i < m_stepElements.size(); i++) {
//     MPointArray aPoints =
//     common::util::swapPointArrayCoordinates(m_stepElements[i-1]->shape,
//     1, 2); MPointArray bPoints =
//     common::util::swapPointArrayCoordinates(m_stepElements[i-0]->shape,
//     1, 2);

//     MPointArray segment;
//     segment.append(aPoints[(int) aPoints.length()/2]);
//     segment.append(bPoints[(int) bPoints.length()/2]);

//     if (common::util::polylinesIntersect(stroke, segment)) {
//       return i;
//     }
//   }
//   return -1;
// }

bool ElementContainer::getClosestPointOnPath(const MPoint pt,
                                             MPoint &closestPt) {
  if (pathLongEnough()) {
    closestPt = common::util::getPointAlongCurve(pt, pathShape.transform());
    return true;
  }
  return false;
}

MPoint ElementContainer::getPointAlongPath(double t) {
  MPoint localPoint = common::util::getPointAlongCurve(t, pathShape.node());

  MFnTransform fnTransform(canvasShape.transform());

  std::ostringstream ss;
  ss << "cmds.xform('";
  ss << fnTransform.fullPathName();
  ss << "', q=True, m=True, ws=True)";
  MString m = MGlobal::executePythonCommandStringResult(ss.str().c_str());
  m = m.substring(1, m.length() - 1);

  MStringArray strVals;
  m.split(',', strVals);
  MDoubleArray vals;
  for (MString s : strVals) {
    vals.append(s.asDouble());
  }

  double dArr[4][4];
  int index = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      dArr[i][j] = vals[index++];
    }
  }
  MMatrix worldXForm(dArr);
  MMatrix canvasMatrix = worldXForm;

  return localPoint * canvasMatrix;
}

MPointArray ElementContainer::getGhostShape() {
  if (!ghostShape.isValid())
    return MPointArray();

  MFnNurbsCurve fnCurve(ghostShape.node());
  MPointArray points;
  fnCurve.getCVs(points);
  return points;
}

} // namespace interface
} // namespace canvas
} // namespace bookmark_squidget