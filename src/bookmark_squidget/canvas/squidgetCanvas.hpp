#ifndef SQUIDGET_CANVAS_HPP
#define SQUIDGET_CANVAS_HPP

#include <memory>
#include <string>
#include <vector>

#include <maya/M3dView.h>
#include <maya/MDagPath.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MSelectionList.h>

#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnTransform.h>

#include "bookmark_squidget/canvas/interface/element.hpp"
#include "bookmark_squidget/canvas/squidget.hpp"
#include "common/maya/factory.hpp"
#include "common/util/math.hpp"
#include "common/util/maya.hpp"

namespace bookmark_squidget {
namespace canvas {
class SquidgetCanvas {
public:
  SquidgetCanvas();
  SquidgetCanvas(MDagPath canvasPath);
  ~SquidgetCanvas();

  void setRestAttributeConfig();
  void deactivate();
  void activate();
  int restAttrSize();

  void HandleCreateModeStroke(MPointArray stroke, MSelectionList selList,
                              bool &shouldDelete);
  bool HandleEditModeStroke(const MPointArray screenStroke, MPoint &pathPoint);

  bool getIntersectingWorldPoint(const MPoint screenPt, MPoint &worldPt);

  void ReadPoint(MPoint point);
  void unhandleSquidget();

  std::vector<std::shared_ptr<interface::StepElement>> m_stepElements;
  std::vector<std::shared_ptr<AttrConfigMap>> m_attributeConfigMaps;
  std::vector<std::shared_ptr<Squidget>> m_squidgets;

  bool isValid();
  void setVisibility(bool isVisible);
  void highlight(bool highlight);
  bool isVisible();
  std::string id;

  std::shared_ptr<AttrConfigMap> rest_attrConfigMap;

private:
  void addSquidget(MIntArray crossingIndices);
  void addStepElement(MPointArray stroke);
  void removeStepElement(int index);
  void addAttributeConfig();
  std::shared_ptr<AttrConfigMap> getAttributeConfig(MStringArray plugs);

  void getSelectedAttributePlugs(MStringArray &arr);
  bool getIntersectingWorldPoints(const MPointArray screenPts,
                                  MPointArray &worldPts);
  MIntArray checkPathCrossing(MPointArray localStroke);
  MIntArray checkElementCrossing(MPointArray stroke);
  bool checkSelfXing(MPointArray localStroke);

  // int getActiveSquidgetIndex();
  int activeSquidgetIndex(MPointArray localPts);
  int _getClosestSquidgetIndex(MPointArray localPts);

  MPointArray localizePointsToCanvas(MPointArray pts);

  std::shared_ptr<interface::StepElement> createStepElement(MPointArray pts);

  MDagPath _canvasPath;

  MSelectionList selectionList;
  MSelectionList restSelectionList;

  bool isHovered;
  int _activeSquidgetIndex = -1;
};
} // namespace canvas
} // namespace bookmark_squidget
#endif // SQUIDGET_CANVAS_HPP