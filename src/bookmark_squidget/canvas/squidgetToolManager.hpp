#ifndef SQUIDGET_TOOL_MANAGER_HPP
#define SQUIDGET_TOOL_MANAGER_HPP

#include <memory>
#include <ostream>
#include <vector>

#include <maya/M3dView.h>
#include <maya/MDagPath.h>
#include <maya/MEulerRotation.h>
#include <maya/MPoint.h>
#include <maya/MPxContext.h>
#include <maya/MSelectionList.h>
#include <maya/MVector.h>

#include <maya/MFnCamera.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnTransform.h>

#include "bookmark_squidget/apps/bookmarkSquidgetTool.hpp"
#include "bookmark_squidget/canvas/squidget.hpp"
#include "bookmark_squidget/canvas/squidgetCanvas.hpp"
#include "common/logging/logger.hpp"

namespace bookmark_squidget {
namespace canvas {
enum SQUIDGET_MODE {
  CREATE,
  EDIT,
  GLOBAL,
};

class SquidgetToolManager {
public:
  SquidgetToolManager();
  ~SquidgetToolManager();

  void createSquidgetCanvas();
  bool toggleSquidgetCanvasesVisibility();
  void updateSquidgetCanvasVisibility(bool visible);

  bool handlePoint(MPoint screenPoint);
  int handlePoints(const MPointArray screenStroke, MPoint &pathPoint);
  void unHandleCanvas();

  int hoverCanvas(short curr_x, short curr_y);
  void canvasHighlight(int index, bool highlight);

  void createMode();
  void editMode();
  void globalMode();
  bool isCreateMode();
  SQUIDGET_MODE mode();

  std::vector<std::shared_ptr<SquidgetCanvas>> _canvases;

private:
  MDagPath createCameraCanvas();

  bool _handlePointsCreate(MPointArray screenStroke);
  int _handlePointsEdit(MPointArray screenStroke, MPoint &pathPoint);
  bool _handlePointEdit(MPoint screenPt);

  int _getClosestSquidgetCanvasIndex(MPoint screenPt);

  void validateTool();

  MSelectionList selList;
  SQUIDGET_MODE m_mode;

  double _canvas_distance;
  double _canvas_width;

  int _hoverCanvasIndex;
  int _currentCanvasIndex;

  bool canvasVisible;

  MString groupName;
};
} // namespace canvas
} // namespace bookmark_squidget

#endif