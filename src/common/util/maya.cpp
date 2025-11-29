#include "common/util/maya.hpp"

namespace common {
namespace util {

MObject getObjectByName(const MString s) {
  MGlobal::selectByName(s, MGlobal::kReplaceList);
  MSelectionList sel_list;
  MGlobal::getActiveSelectionList(sel_list);
  MDagPath path;
  sel_list.getDagPath(0, path);
  MObject outputTrans = path.transform();
  return outputTrans;
}

MObject resampleCurve(MObject curve, int numSamples) {
  MFnDagNode fnDagNode(curve);
  MString curveName = fnDagNode.partialPathName();
  std::ostringstream ss;
  ss << "rebuildCurve -ch 0 -rpo 1 -rt 0 -end 1 -kr 0 -kcp 0 -kep 1 -kt 0 -s";
  ss << " " << numSamples;
  ss << " -d 3 -tol 0.01 ";
  ss << "\"" << curveName << "\"";
  std::string cmd = ss.str();
  // cout << "command: " << cmd << endl;

  MGlobal::executeCommand(MString(cmd.c_str())); // replaces original

  return curve;
  // IDK How this works.  Documentation is bad.
  // MObject newCurveData = fnCurve.rebuild(numSamples, DEGREE,
  //                                     KEEP_RANGE, END_KNOTS,
  //                                     KEEP_END_POINTS, KEEP_TANGENTS,
  //                                     KEEP_CONTROL_POINTS, TOLERANCE,
  //                                     &status);
}

MStatus getObjectAttributePair(const MString s, MDagPath &dagPath,
                               MPlug &plug) {
  MStringArray chArr;
  s.split('.', chArr);
  MString obj = chArr[0];
  MString attr = chArr[1];

  // Save Prev Selection
  MSelectionList prevSel;
  MGlobal::getActiveSelectionList(prevSel);

  // Get selection obj path
  MSelectionList sel;
  MGlobal::selectByName(obj, MGlobal::kReplaceList);
  MGlobal::getActiveSelectionList(sel);

  // Get plug
  sel.getDagPath(0, dagPath);
  MFnDagNode fnDagNode;
  fnDagNode.setObject(dagPath);
  MPlugArray plugs;
  fnDagNode.getConnections(plugs);
  plug = fnDagNode.findPlug(attr, true);

  MGlobal::setActiveSelectionList(prevSel);
  return MS::kSuccess;
}

MPointArray worldToScreen(const MPointArray worldPoints) {
  M3dView view = M3dView::active3dView();
  MPointArray screenPoints;
  for (MPoint p : worldPoints) {
    MPoint screenPoint;
    short x, y;
    view.worldToView(p, x, y);
    screenPoints.append(MPoint(x, y, 0));
  }
  return screenPoints;
}

Eigen::Matrix4Xd pointArrayToMatrix(const MPointArray &points) {
  Eigen::Matrix4Xd ret(4, points.length());
  for (unsigned int i = 0; i < points.length(); ++i) {
    ret.col(i) << points[i].x, points[i].y, points[i].z, 1;
  }
  return ret;
}

MPointArray MatrixToPointArray(const Eigen::Matrix4Xd &points) {
  MPointArray ret;
  for (int i = 0; i < points.cols(); ++i) {
    ret.append(MPoint(points(0, i), points(1, i), points(2, i), 1));
  }
  return ret;
}

Eigen::Matrix4d mMatrixToMatrix(const MMatrix &m) {
  Eigen::Matrix4d ret;
  for (int i = 0; i < 4; ++i) {
    ret.row(i) << m[i][0], m[i][1], m[i][2], m[i][3];
  }
  return ret;
}

MMatrix matrixToMMatrix(const Eigen::Matrix4d &m) {
  MMatrix ret;
  for (int i = 0; i < 4; ++i) {
    ret[i][0] = m(i, 0);
    ret[i][1] = m(i, 1);
    ret[i][2] = m(i, 2);
    ret[i][3] = m(i, 3);
  }
  return ret;
}

MPointArray swapPointArrayCoordinates(MPointArray arr, int a, int b) {
  MPointArray points;
  if ((a == 1 && b == 2) || (a == 2 && b == 1)) {
    for (MPoint p : arr) {
      points.append(MPoint(p.x, p.z, p.y));
    }
  }

  return points;
}

MObjectArray combineAndFilterCurves(MObjectArray curves) {
  MObjectArray ret;
  // Keep visible points, combine curve ends

  for (MObject curve : curves) {
    MFnNurbsCurve fnCurve(curve);
    MPointArray points;
    fnCurve.getCVs(points, MSpace::kWorld);

    MPointArray visiblePoints = filterInvisiblePoints(points);
    if (visiblePoints.length() == 0) {
      continue;
    }
    fnCurve.setCVs(visiblePoints, MSpace::kWorld);
  }

  return ret;
}

MPointArray filterInvisiblePoints(MPointArray points) {
  MPointArray ret;

  for (MPoint p : points) {
    if (!pointIntersection(p)) {
      ret.append(p);
    }
  }
  return ret;
}

bool pointIntersection(MPoint pointToCheck) {
  // Get the camera position and direction
  M3dView view = M3dView::active3dView();
  short x, y;
  view.worldToView(pointToCheck, x, y);
  MPoint cameraPosition;
  MVector cameraDirection;
  view.viewToWorld(x, y, cameraPosition, cameraDirection);
  cameraDirection = MVector(pointToCheck - cameraPosition);

  // Iterate through all meshes in the scene
  bool pointOccluded = false;
  MItDag dagIterator(MItDag::kDepthFirst, MFn::kMesh);
  for (; !dagIterator.isDone(); dagIterator.next()) {
    MDagPath path;
    dagIterator.getPath(path);
    MFnMesh mesh(path);

    MFloatPoint hitPoint;
    float param;
    // F this method.  Check return status.
    bool intersects = mesh.closestIntersection(
        cameraPosition, cameraDirection, nullptr, nullptr, false,
        MSpace::kWorld, 10000, false, nullptr, hitPoint, &param, nullptr,
        nullptr, nullptr, nullptr, 1e-6);

    // Check for ray intersection with the mesh
    if (intersects && param < .99) {
      pointOccluded |= intersects;
    }
  }
  return pointOccluded;
}

void deleteCurve(MObject curve) {
  MDagPath path;
  MFnDagNode fnDagNode(curve);
  fnDagNode.getPath(path);
  MObject curveTrans = path.transform();
  MGlobal::deleteNode(curveTrans);
}

MString convertStringArrayToString(const MStringArray &stringArray,
                                   const MString &delimiter) {
  MString result;

  for (unsigned int i = 0; i < stringArray.length(); ++i) {
    result += stringArray[i];

    // Add the delimiter if it's not the last element
    if (i < stringArray.length() - 1) {
      result += delimiter;
    }
  }

  return result;
}

MDagPath getDagPathFromName(const MString &name) {
  MObject geomObj = getObjectByName(name);
  MFnDagNode geomFn(geomObj);
  MDagPath geomPath;
  geomFn.getPath(geomPath);
  return geomPath;
}

MDagPath getNurbsCurveDagPathFromObject(MObject &object) {
  MStatus status;
  MFnDagNode fnDagNode(object);
  MDagPath curvePath;
  status = fnDagNode.getPath(curvePath);
  CHECK_MSTATUS(status);
  return curvePath;
}

MPoint getWorldToScreenPoint(const MPoint &worldPt) {
  M3dView view = M3dView::active3dView();
  short x, y;
  view.worldToView(worldPt, x, y);
  MPoint mayaP = MPoint(x, y, 0);

  return mayaP;
}
} // namespace util
} // namespace common