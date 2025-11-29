#include "common/util/math.hpp"

namespace common {
namespace util {

int polylinesIntersect(MPointArray polyline1, MPointArray polyline2) {
  // Check if polylines intersect
  int times = 0;
  for (size_t i = 0; i < polyline1.length() - 1; i++) {
    for (size_t j = 0; j < polyline2.length() - 1; j++) {
      Pt p1 = polyline1[i];
      Pt p2 = polyline1[i + 1];
      Pt p3 = polyline2[j];
      Pt p4 = polyline2[j + 1];

      if (intersect(p1, p2, p3, p4)) {
        times++;
      }
    }
  }

  return times;
}

bool intersect(Pt A, Pt B, Pt C, Pt D) {
  return ccw(A, C, D) != ccw(B, C, D) && ccw(A, B, C) != ccw(A, B, D);
}

bool ccw(Pt A, Pt B, Pt C) {
  return (C.y - A.y) * (B.x - A.x) > (B.y - A.y) * (C.x - A.x);
}

MPoint getPointAlongCurve(MPoint pt, MObject curve) {
  double param = getDistanceAlongCurve(pt, curve);
  return getPointAlongCurve(param, curve);
}

MPoint getPointAlongCurve(double param, MObject curve) {
  MPoint pt;

  MFnNurbsCurve curveFn(curve);
  MPointArray pts;
  curveFn.getCVs(pts);
  curveFn.getPointAtParam(param, pt);

  return pt;
}

double getDistanceAlongCurve(MPoint pt, MObject curve) {
  double param;
  MFnNurbsCurve curveFn(curve);
  MPoint closestPoint =
      curveFn.closestPoint(pt, &param, 0.0001, MSpace::kObject);
  return param;
}

MPointArray shiftPoints(MPointArray points, MPoint newOrigin) {
  MPointArray newPoints;
  for (MPoint pt : points) {
    MPoint newPt = pt - newOrigin;
    newPoints.append(newPt);
  }
  return newPoints;
}

Eigen::VectorXd nurbsToVector(MObject curve) {
  MFnNurbsCurve fnNurbs(curve);

  MPointArray cvs;
  fnNurbs.getCVs(cvs);
  // cout << cvs << endl;

  Eigen::VectorXd vec(3 * cvs.length());
  for (unsigned int i = 0; i < cvs.length(); ++i) {
    MPoint p = cvs[i];
    vec[3 * i + 0] = p.x;
    vec[3 * i + 1] = p.y;
    vec[3 * i + 2] = p.z;
  }
  return vec;
}

std::chrono::milliseconds getTime() {
  using namespace std::chrono;
  milliseconds ms =
      duration_cast<milliseconds>(system_clock::now().time_since_epoch());
  return ms;
}

std::string generateRandomID(int length) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> distribution(0, 9);

  // Generate random digits
  std::ostringstream oss;
  for (int i = 0; i < length; ++i) {
    oss << distribution(gen);
  }

  return oss.str();
}

} // namespace util
} // namespace common