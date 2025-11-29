#ifndef ELEMENT_CONTAINER_HPP
#define ELEMENT_CONTAINER_HPP

#include <memory>
#include <vector>

#include <Eigen/Dense>

#include <maya/MDagPath.h>
#include <maya/MGlobal.h>
#include <maya/MObjectArray.h>
#include <maya/MPlug.h>
#include <maya/MPointArray.h>
#include <maya/MSelectionList.h>
#include <maya/MVector.h>

#include <maya/MFnBlendShapeDeformer.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnSphereData.h>
#include <maya/MFnTransform.h>

#include "bookmark_squidget/canvas/interface/element.hpp"
#include "common/maya/factory.hpp"
#include "common/util/math.hpp"
#include "common/util/maya.hpp"

namespace bookmark_squidget {
namespace canvas {
namespace interface {

/**
 * @brief Handles all physical elements that compose the squidget.
 * Will add and remove step curves and queries t value along the path.
 */
class ElementContainer {
public:
  ElementContainer();
  ElementContainer(MDagPath canvasPath);
  ~ElementContainer();

  void AddStepElement(std::shared_ptr<StepElement> stepElement);
  // void AddStepElement(MPointArray stroke);
  // void DeleteStepElement(int index);

  double UpdateValueAlongPath(MPoint pt);
  int CheckStrokeCrossing(MPointArray stroke);
  double CheckPathCrossing(MPointArray stroke);

  double queryStrokeValue(const MPointArray localPts, double &dist);

  bool setTValue(double t);
  MPoint getPointAlongPath(double t);

  std::vector<std::shared_ptr<StepElement>> m_stepElements;
  std::vector<double> m_stepElementKey;

  MDagPath pathShape;
  MDagPath ghostShape;
  MDagPath canvasShape;

  MFnBlendShapeDeformer deform;
  MObject squidget_t;
  MPointArray getGhostShape();

private:
  void createInterface();

  void createPath();
  void createTAttribute();
  void createGhost();
  void createElementBlendShape();
  void connectSquidgetToBlendshape();

  void updateInterface();
  void deleteInterface();

  bool pathLongEnough();
  MPointArray localizePointsToCanvas(MPointArray pts);
  std::shared_ptr<StepElement> createStepElement(MPointArray pts);

  bool getClosestPointOnPath(const MPoint pt, MPoint &closestPt);

  Eigen::VectorXd translateToOrigin(Eigen::VectorXd f);
  double curveDistance(Eigen::VectorXd a, Eigen::VectorXd b);
};

} // namespace interface
} // namespace canvas
} // namespace bookmark_squidget

#endif