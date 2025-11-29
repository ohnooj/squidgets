#ifndef UTIL_MAYA_HPP
#define UTIL_MAYA_HPP

#include <Eigen/Dense>
#include <sstream>

#include <maya/M3dView.h>
#include <maya/MDagPath.h>
#include <maya/MDataHandle.h>
#include <maya/MGlobal.h>
#include <maya/MMatrix.h>
#include <maya/MObjectArray.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MPointArray.h>
#include <maya/MSelectionList.h>
#include <maya/MString.h>

#include <maya/MFnMesh.h>
#include <maya/MItDag.h>

#include <maya/MFnDagNode.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsCurveData.h>

namespace common {
namespace util {

#define DEGREE 1
#define KEEP_RANGE 0
#define END_KNOTS 1
#define KEEP_END_POINTS 1
#define KEEP_TANGENTS 0
#define KEEP_CONTROL_POINTS 0
#define TOLERANCE 0.01

MObject resampleCurve(MObject curve, int numSamples);
MPointArray worldToScreen(const MPointArray worldPoints);

Eigen::Matrix4Xd pointArrayToMatrix(const MPointArray &points);
MPointArray MatrixToPointArray(const Eigen::Matrix4Xd &points);

Eigen::Matrix4d mMatrixToMatrix(const MMatrix &m);
MMatrix matrixToMMatrix(const Eigen::Matrix4d &m);

MPointArray swapPointArrayCoordinates(MPointArray arr, int a, int b);

MObjectArray combineAndFilterCurves(MObjectArray curves);
MPointArray filterInvisiblePoints(MPointArray points);
bool pointIntersection(MPoint pointToCheck);

void deleteCurve(MObject curve);
MString convertStringArrayToString(const MStringArray &stringArray,
                                   const MString &delimiter = "");

// Getters
MObject getObjectByName(const MString s);
MDagPath getDagPathFromName(const MString &name);
MDagPath getNurbsCurveDagPathFromObject(MObject &object);
MStatus getObjectAttributePair(const MString s, MDagPath &dagPath, MPlug &plug);
MPoint getWorldToScreenPoint(const MPoint &worldPt);
} // namespace util
} // namespace common

#endif
