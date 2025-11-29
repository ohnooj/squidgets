#include "bookmark_squidget/canvas/squidgetToolManager.hpp"

namespace bookmark_squidget {
namespace canvas {
SquidgetToolManager::SquidgetToolManager() {
  // cout << "SquidgetToolManager()" << std::endl;
  m_mode = SQUIDGET_MODE::CREATE;
  _canvas_distance = 13;
  _canvas_width = 7;
  _currentCanvasIndex = -1;
  _hoverCanvasIndex = -1;
  canvasVisible = true;

  // Clear SQUIDGETS group or create if needed
  groupName = MString("SQUIDGETS");
  std::ostringstream ss;
  ss << "squidgets_group = cmds.ls('" << groupName << "')" << endl;
  ss << "if len(squidgets_group) == 0:" << endl;
  ss << "    cmds.group(em=True, n='" << groupName << "')" << endl;
  ss << "else:" << endl;
  ss << "    group_children = cmds.listRelatives('" << groupName << "') or []"
     << endl;
  MGlobal::executePythonCommand(ss.str().c_str());

  // Create canvas color material
  ss.str("");
  ss << "def canvasMat(name, t):" << endl;
  ss << "    canvasHL = cmds.shadingNode('lambert', asShader=True, n=name)"
     << endl;
  ss << "    canvasHLSG = f'{canvasHL}SG'" << endl;
  ss << "    cmds.sets(canvasHL, n=canvasHLSG, renderable=True, "
        "noSurfaceShader=True, empty=True)"
     << endl;
  ss << "    cmds.connectAttr(f'{canvasHL}.outColor', "
        "f'{canvasHLSG}.surfaceShader', f=True)"
     << endl;
  ss << "    cmds.setAttr(f'{canvasHL}.transparency', t, t, t, type='double3') "
     << endl;
  ss << "    cmds.setAttr(f'{canvasHL}.color', 0.75, 0.75, 0.0, type='double3')"
     << endl;
  ss << "    cmds.setAttr(f'{canvasHL}.ambientColor', 1.0, 1.0, 1.0, "
        "type='double3')"
     << endl;
  ss << "    return canvasHL" << endl;
  MGlobal::executePythonCommand(ss.str().c_str());

  ss.str("");
  ss << "canvasNormalExists = cmds.ls('canvasNormalColor')" << endl;
  ss << "if len(canvasNormalExists) == 0:" << endl;
  ss << "    canvasMat('canvasNormalColor', 0.94)" << endl;
  ss << "    canvasMat('canvasHighlightColor', 0.9)" << endl;
  MGlobal::executePythonCommand(ss.str().c_str());
}

SquidgetToolManager::~SquidgetToolManager() {
  selList.clear();
}

int SquidgetToolManager::handlePoints(const MPointArray screenStroke,
                                      MPoint &pathPoint) {
  /**
   * Input: 2D curvePoints
   */

  MGlobal::getActiveSelectionList(selList);

  // cout << "SquidgetToolManager::handlePoints()" << endl;
  if (mode() == SQUIDGET_MODE::CREATE) {
    return _handlePointsCreate(screenStroke);
  }
  if (mode() == SQUIDGET_MODE::EDIT) {
    return _handlePointsEdit(screenStroke, pathPoint);
  }
  return 1; // On border
}

bool SquidgetToolManager::handlePoint(MPoint screenPoint) {
  if (mode() == SQUIDGET_MODE::EDIT) {
    return _handlePointEdit(screenPoint);
  }
  return false;
}

void SquidgetToolManager::unHandleCanvas() {
  if (_currentCanvasIndex >= 0) {
    std::string canvID = _canvases.at(_currentCanvasIndex)->id;
  }
  for (auto canvas : _canvases) {
    canvas->unhandleSquidget();
  }

  _currentCanvasIndex = -1;
}

void SquidgetToolManager::validateTool() {
  /**
   * @brief Validate the tool if things were deleted.
   */

  unsigned int i = 0;
  while (i < _canvases.size()) {
    if (!_canvases[i]->isValid()) {
      _canvases.erase(_canvases.begin() + i);
    } else {
      ++i;
    }
  }
}

// =============================================================================
// CREATION MODE
// =============================================================================
bool SquidgetToolManager::_handlePointsCreate(MPointArray screenStroke) {
  std::shared_ptr<SquidgetCanvas> activeCanvas;
  int index = _getClosestSquidgetCanvasIndex(screenStroke[0]);
  if (index == -1) {
    return false;
  }
  activeCanvas = _canvases.at(index);

  bool shouldDelete = false;
  activeCanvas->HandleCreateModeStroke(screenStroke, selList, shouldDelete);

  if (shouldDelete) {
    // Find activeCanvas in _canvases and delete it
    auto it = std::find(_canvases.begin(), _canvases.end(), activeCanvas);
    _canvases.erase(it);
    _currentCanvasIndex = -1;
    _hoverCanvasIndex = -1;
  }

  return true;
}

int SquidgetToolManager::_getClosestSquidgetCanvasIndex(MPoint screenPt) {
  if (_canvases.size() < 1) {
    return -1;
  }

  MPoint worldPt;
  MVector worldVec;
  M3dView view = M3dView::active3dView();
  view.viewToWorld(screenPt.x, screenPt.y, worldPt, worldVec);

  double minDist = std::numeric_limits<double>::max();
  int idx = -1;
  for (size_t i = 0; i < _canvases.size(); ++i) {
    auto can = _canvases[i];

    // skip invisibile canvases
    if (!can->isVisible()) {
      continue;
    }

    MPoint intersectPt;
    bool intersect = can->getIntersectingWorldPoint(screenPt, intersectPt);
    if (intersect) {
      double dist = worldPt.distanceTo(intersectPt);
      if (dist < minDist) {
        minDist = dist;
        idx = i;
      }
    }
  }

  return idx;
}

// Canvas Methods ==============================================================
void SquidgetToolManager::createSquidgetCanvas() {
  MSelectionList selList;
  MGlobal::getActiveSelectionList(selList);

  if (selList.length() < 1) {
    return;
  }

  MDagPath canvasPath = createCameraCanvas();

  MObject canvasTransform = canvasPath.transform();
  MFnDagNode fnDagNode;
  fnDagNode.setObject(canvasTransform);

  std::ostringstream ss;
  ss << "parent ";
  ss << fnDagNode.fullPathName() << " ";
  ss << "SQUIDGETS;";
  MGlobal::executeCommand(ss.str().c_str());

  MDagPath freeCanvasPath;
  MDagPath canvasTransPath;
  fnDagNode.setObject(canvasTransform);
  fnDagNode.getPath(freeCanvasPath);
  fnDagNode.getPath(canvasTransPath);
  freeCanvasPath.extendToShape();

  MGlobal::setActiveSelectionList(selList);
  std::shared_ptr<SquidgetCanvas> activeCanvas;
  activeCanvas = std::make_shared<SquidgetCanvas>(freeCanvasPath);
  if (activeCanvas->restAttrSize() == 0) {
    ss.str("");
    ss.clear();
    ss << "cmds.delete('" << canvasTransPath.fullPathName() << "')" << endl;
    MGlobal::executePythonCommand(ss.str().c_str());
    return;
  }

  _canvases.push_back(activeCanvas);

  MStringArray selNames;
  selList.getSelectionStrings(selNames);
  MString allSelNames;
  for (unsigned int i = 0; i < selNames.length(); ++i) {
    allSelNames += selNames[i];
    if (i < selNames.length() - 1) {
      allSelNames += ", ";
    }
  }

  MString canvasName = canvasTransPath.partialPathName() + ": " + allSelNames;
  ss.str("");
  ss.clear();
  ss << "annot = cmds.annotate('" << canvasTransPath.fullPathName()
     << "', tx=' " << canvasName << "');" << endl;
  ss << "cmds.setAttr(annot + '.displayArrow', 0);" << endl;
  ss << "annot = cmds.listRelatives(annot, p=True, f=True)[0]" << endl;
  ss << "cmds.move(0, 0, " << -_canvas_width / 2 << ", annot)" << endl;
  ss << "cmds.parent(annot, '" << canvasTransPath.fullPathName() << "', r=True)"
     << endl;
  MGlobal::executePythonCommand(ss.str().c_str());

  MGlobal::clearSelectionList();
}

int SquidgetToolManager::hoverCanvas(short curr_x, short curr_y) {
  // cout << "SquidgetToolmanager::HoverCanvas" << endl;
  int currHoverCanvas = _getClosestSquidgetCanvasIndex(MPoint(curr_x, curr_y));
  // cout << "Curr Hover Canvas: " << currHoverCanvas << endl;
  // cout << "HoverCanvas canvas: " << _hoverCanvasIndex << endl;

  if (_hoverCanvasIndex != currHoverCanvas) {
    canvasHighlight(currHoverCanvas, true);
    // cout << "Highlight: " << currHoverCanvas << endl;
    canvasHighlight(_hoverCanvasIndex, false);
    // cout << "UnHighlight: " << _hoverCanvasIndex<< endl;
  }
  // cout << "after" << endl;
  _hoverCanvasIndex = currHoverCanvas;
  return _hoverCanvasIndex;
}

void SquidgetToolManager::canvasHighlight(int index, bool highlight) {
  if (index < 0) {
    return;
  }

  std::shared_ptr<SquidgetCanvas> activeCanvas = _canvases.at(index);
  activeCanvas->highlight(highlight);
}

MDagPath SquidgetToolManager::createCameraCanvas() {
  /**
   * @brief Creates a canvas and attaches it to the current active camera
   * viewport.
   */
  double w = _canvas_width;

  // Create nurbs plane
  std::ostringstream ss;
  ss << "$cn = `nurbsPlane -p 0 0 0 -ax 0 1 0 -w " << w
     << " -lr 1 -d 3 -u 1 -v 1 -ch 1 ";
  ss << "-n \"canvas\" -o true`;";
  ss << "$res = $cn[0];";
  MString res = MGlobal::executeCommandStringResult(ss.str().c_str());
  // cout << "Command: " << ss.str() << endl;

  MStringArray nurbsPlaneRes;
  res.split(';', nurbsPlaneRes);
  MString canvasName = nurbsPlaneRes[0];
  MString canvasCreateNode = nurbsPlaneRes[1];

  ss.str("");
  ss << "print('" << canvasName << "')" << endl;
  ss << "cmds.sets('" << canvasName
     << "', e=True, forceElement=\"canvasNormalColorSG\");";
  // cout << "Command: " << ss.str() << endl;
  MGlobal::executePythonCommand(ss.str().c_str());

  // Set current _canvasPath to nurbs plane
  MGlobal::selectByName(canvasName, MGlobal::kReplaceList);
  MSelectionList selection;
  MGlobal::getActiveSelectionList(selection);

  MDagPath canvasPath;
  selection.getDagPath(0, canvasPath);
  MObject canvasTransform = canvasPath.transform();

  // Set canvas as child of active camera to move with camera
  MDagPath cameraPath;
  M3dView view = M3dView::active3dView();
  view.getCamera(cameraPath);
  MObject cameraTransform = cameraPath.transform();

  MFnDagNode fnDagNode;
  fnDagNode.setObject(cameraTransform);
  fnDagNode.addChild(canvasTransform);

  // Transform plane to face viewport
  MFnCamera fnCamera(cameraPath);
  MPoint camEye = fnCamera.eyePoint();
  MVector camDir = fnCamera.viewDirection(MSpace::kWorld);

  MFnTransform fnTransform(canvasPath);
  fnTransform.translateBy(_canvas_distance * camDir, MSpace::kWorld);
  fnTransform.rotateBy(MEulerRotation(M_PI / 2, 0, 0));
  // cout << "What is res: " << res << endl;

  ss.str("");
  ss.clear();
  ss << "print('" << canvasPath.fullPathName() << "')" << endl;
  ss << "cmds.setAttr('" << canvasPath.fullPathName() << ".overrideEnabled', 1)"
     << endl;
  ss << "cmds.setAttr('" << canvasPath.fullPathName()
     << ".overrideDisplayType', 2)" << endl;
  MGlobal::executePythonCommand(ss.str().c_str());
  // cout << "Command: " << ss.str() << endl;
  return canvasPath;
}

// =============================================================================
// EDIT MODE
// =============================================================================
int SquidgetToolManager::_handlePointsEdit(MPointArray screenStroke,
                                           MPoint &pathPoint) {

  bool hasCanvas = _currentCanvasIndex >= 0;
  bool allSameIndex = false;
  // Use stroke to check a canvas if no canvas is selected.
  if (!hasCanvas) {
    int index = _getClosestSquidgetCanvasIndex(screenStroke[0]);

    allSameIndex =
        std::all_of(screenStroke.begin(), screenStroke.end(), [&](MPoint pt) {
          return _getClosestSquidgetCanvasIndex(pt) == index;
        });

    _currentCanvasIndex = (index == -1 || !allSameIndex) ? -1 : index;
  }

  hasCanvas = _currentCanvasIndex >= 0;

  if (hasCanvas) {
    std::shared_ptr<SquidgetCanvas> activeCanvas =
        _canvases.at(_currentCanvasIndex);
    return activeCanvas->HandleEditModeStroke(screenStroke, pathPoint);
  } else if (allSameIndex) {
    return 2;
  } else {
    // Crossing canvas border
    // Stroke must be on canvas completely: Squidget Edit
    // Or off completely: Global Mode
    return 1;
  }
}

bool SquidgetToolManager::_handlePointEdit(MPoint screenPt) {
  MGlobal::getActiveSelectionList(selList);
  if (mode() == SQUIDGET_MODE::EDIT) {
    // Assume canvas index is valid
    if (_currentCanvasIndex >= 0) {
      _canvases.at(_currentCanvasIndex)->ReadPoint(screenPt);
      return true;
    }
  }
  return false;
}

bool SquidgetToolManager::toggleSquidgetCanvasesVisibility() {
  if (mode() == SQUIDGET_MODE::CREATE) {
    return false;
  }

  if (_canvases.size() < 1) {
    return true;
  }

  canvasVisible = !canvasVisible;
  updateSquidgetCanvasVisibility(canvasVisible);
  return canvasVisible;
}

void SquidgetToolManager::updateSquidgetCanvasVisibility(bool visible) {
  for (auto canvas : _canvases) {
    canvas->setVisibility(visible);
  }
}

// =============================================================================
// TOGGLE MODE
// =============================================================================
void SquidgetToolManager::createMode() {
  m_mode = SQUIDGET_MODE::CREATE;
  unHandleCanvas();
  for (auto canvas : _canvases) {
    canvas->deactivate();
  }

  canvasVisible = true;
  updateSquidgetCanvasVisibility(canvasVisible);
}

void SquidgetToolManager::editMode() {
  m_mode = SQUIDGET_MODE::EDIT;
  unHandleCanvas();

  for (auto canvas : _canvases) {
    canvas->activate();
  }

  canvasVisible = true;
  updateSquidgetCanvasVisibility(canvasVisible);
}

SQUIDGET_MODE SquidgetToolManager::mode() { return m_mode; }

bool SquidgetToolManager::isCreateMode() {
  return m_mode == SQUIDGET_MODE::CREATE;
}

} // namespace canvas
} // namespace bookmark_squidget