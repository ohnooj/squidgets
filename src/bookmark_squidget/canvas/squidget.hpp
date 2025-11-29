#ifndef SQUIDGET_HPP
#define SQUIDGET_HPP

#include <iostream>

#include "maya/M3dView.h"
#include "maya/MDagPath.h"
#include "maya/MDoubleArray.h"
#include "maya/MGlobal.h"
#include "maya/MPointArray.h"
#include "maya/MString.h"

#include "maya/MFnNurbsCurve.h"
#include "maya/MFnNurbsSurface.h"

#include "bookmark_squidget/canvas/attributeConfigManager.hpp"
#include "bookmark_squidget/canvas/interface/elementContainer.hpp"
#include "common/logging/logger.hpp"
#include "common/maya/factory.hpp"
#include "common/util/math.hpp"

namespace bookmark_squidget {
namespace canvas {
class Squidget {
  /**
   * @brief How squidget works.
   * Squidet takes in a stroke.  We determin the t value [0, 1] of the stroke.
   * We then use the t value to drive two things:
   * - Element Shape: The element shapes will have a blendshape node with
   * inbetween shapes.
   * - Attribute Configuration interpolation:
   *     - AttributeConfiguration will have a special attr that uses
   * setDrivenKey with keys [config1, config2, config3]
   */
public:
  Squidget();
  Squidget(MDagPath canvasPath);
  ~Squidget();

  void addStepElement(std::shared_ptr<interface::StepElement> stepElement);
  void addAttributeConfig(std::shared_ptr<AttrConfigMap> attrConfigMap);

  void activate();
  void deactivate();

  void SetValueFromStroke(const MPointArray screenStroke, MPoint &pathPoint);
  double getTValue() { return t; }
  bool setTValue(double t);

  void ReadPoint(MPoint point);

  bool checkStrokeCrossing(MPointArray localPts);
  double estimateStrokeValue(const MPointArray worldPts, double &dist);
  void setSingleStrokeAttributes();

  MObject getPathShape();
  void setSingleStroke(bool singleStroke) { isSingleStroke = singleStroke; }
  bool getSingleStroke() { return isSingleStroke; }

  std::shared_ptr<canvas::AttributeConfigManager> attributeConfigManager;
  std::shared_ptr<canvas::interface::ElementContainer> elementContainer;
  std::string id;

private:
  void setDrivenConfigs();
  void clearDrivenKeyframes();

  void setDiscrete();

  MSelectionList selectionList;
  MStringArray drivenKeyframes;

  MPoint cam;
  MVector view;
  double t;
  MPointArray lastStroke;
  bool isSingleStroke = false;
};

} // namespace canvas
} // namespace bookmark_squidget
#endif
