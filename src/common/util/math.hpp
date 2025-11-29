#ifndef MATH_HPP
#define MATH_HPP

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <utility> // For std::pair
#include <vector>

#include <Eigen/Dense>

#include <maya/MDagPath.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>

#include <maya/MFnNurbsCurve.h>

namespace common {
namespace util {
namespace {
using Pt = MPoint;
}

int polylinesIntersect(MPointArray polyline1, MPointArray polyline2);
bool intersect(Pt A, Pt B, Pt C, Pt D);
bool ccw(Pt A, Pt B, Pt C);

MPoint getPointAlongCurve(MPoint pt, MObject curve);
MPoint getPointAlongCurve(double param, MObject curve);
double getDistanceAlongCurve(MPoint pt, MObject curve);

MPointArray shiftPoints(MPointArray points, MPoint newOrigin);

Eigen::VectorXd nurbsToVector(MObject curve);

std::chrono::milliseconds getTime();
std::string generateRandomID(int length);

} // namespace util
} // namespace common

#endif