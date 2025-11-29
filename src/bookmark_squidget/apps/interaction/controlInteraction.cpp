#include "bookmark_squidget/apps/interaction/controlInteraction.hpp"

namespace bookmark_squidget {
namespace apps {
namespace interaction {
ControlInteraction::ControlInteraction(
    std::shared_ptr<canvas::SquidgetToolManager> squidgetToolManager) {
  this->squidgetToolManager = squidgetToolManager;
  callbackId = -1;

  setUI();
  editDragMode = false;
  queryMode = false;
  this->squidgetToolManager->editMode();
}

void ControlInteraction::modeSetup() { removeCallback(); }

void ControlInteraction::modeCleanup() { removeCallback(); }

MStatus ControlInteraction::doPress(MEvent &event,
                                    MHWRender::MUIDrawManager &drawMgr,
                                    const MHWRender::MFrameContext &context) {
  // cout << "ControlInteraction::doPress" << endl;
  MGlobal::getActiveSelectionList(prevSelection);
  event.getPosition(curr_x, curr_y);

  screenStrokePoints.clear();
  screenStrokePoints.append(MPoint(curr_x, curr_y, 0));

  dragTime = common::util::getTime();
  last_x = curr_x;
  last_y = curr_y;

  startCallback(curr_x, curr_y, &drawMgr);
  if (event.isModifierShift()) { // Set query stroke mode
    queryMode = true;            // May need to unset somewhere else
  } else {
    queryMode = false;
  }

  return MS::kSuccess;
}

MStatus ControlInteraction::doDrag(MEvent &event,
                                   MHWRender::MUIDrawManager &drawMgr,
                                   const MHWRender::MFrameContext &context) {
  // cout << "ControlInteraction::doDrag" << endl;
  short x, y;
  event.getPosition(x, y);
  return doDrag(x, y, event.isModifierControl(), drawMgr);
}

MStatus ControlInteraction::doDrag(short x, short y, bool isCtrl,
                                   MHWRender::MUIDrawManager &drawMgr) {
  // cout << "ControlInteraction::doDrag helper " << editDragMode << endl;
  curr_x = x;
  curr_y = y;

  double distFromLast = sqrt(pow(curr_x - last_x, 2) + pow(curr_y - last_y, 2));
  if (distFromLast > strokePointDist) {
    screenStrokePoints.append(MPoint(curr_x, curr_y, 0));
    last_x = curr_x;
    last_y = curr_y;
  }

  if (!editDragMode) {
    startCallback(x, y, &drawMgr);
  }

  if (editDragMode) {
    removeCallback();
    MPoint mouse(curr_x, curr_y, 0);
    MPoint offsetMouse(curr_x + offsetVec.x, curr_y + offsetVec.y, offsetVec.z);
    bool usedSquidget = squidgetToolManager->handlePoint(offsetMouse);
    if (!usedSquidget) {
      tool->ReadPoint(mouse, offsetVec, isCtrl);
    }
  }
  return MS::kSuccess;
}

void ControlInteraction::removeCallback() {
  if (callbackId != static_cast<MCallbackId>(0)) {
    MTimerMessage::removeCallback(callbackId);
  }
  callbackId = static_cast<MCallbackId>(0);
}

void ControlInteraction::startCallback(short x, short y,
                                       MHWRender::MUIDrawManager *drawMgr) {
  removeCallback();
  MyTimerCallbackData *callbackData =
      new MyTimerCallbackData(this, x, y, drawMgr);
  callbackId =
      MTimerMessage::addTimerCallback(.4f, timerCallback, callbackData);
}

void ControlInteraction::timerCallback(float elapsedTime, float lastTime,
                                       void *clientData) {
  MyTimerCallbackData *obj = static_cast<MyTimerCallbackData *>(clientData);
  if (obj && obj->obj && obj->uiDrawManager && !obj->obj->isDragMode()) {
    obj->obj->setDragMode();
    // cout << "Is drag mode " << obj->obj->isDragMode() << endl;
    if (!obj->obj->isDragMode()) {
      return;
    }
    // cout << obj->x << " " << obj->y << endl;
    obj->obj->doDrag(obj->x, obj->y, false, *obj->uiDrawManager);
    obj->obj->draw(*(obj->uiDrawManager));

    M3dView view = M3dView::active3dView();
    view.refresh();
    // cout << "elapsed Time: " << elapsedTime << " seconds, lasttime: " <<
    // lastTime << endl;
  }
}

bool ControlInteraction::isDragMode() { return editDragMode; }

void ControlInteraction::setDragMode() {
  MPoint controlPt = MPoint::origin;

  int handlePointsCode =
      squidgetToolManager->handlePoints(screenStrokePoints, controlPt);
  if (handlePointsCode == 2) { // All off canvas.
    tool->setPenControlPoints(screenStrokePoints);
    // BookmarkSquidgetTool.parseArgs()
    tool->queryClosestObjectOnly();
    tool->predoIt();
    tool->redoIt();
    tool->getControlPt(controlPt);
    if (controlPt == MPoint::origin) {
      // cout << "Nothing in setDragMode" << endl;
      return;
    }
  }

  MMatrix cameraMatrix, projMatrix, screenMatrix;
  M3dView view = M3dView::active3dView();
  int width = view.portWidth();
  int height = view.portHeight();
  view.modelViewMatrix(cameraMatrix);
  view.projectionMatrix(projMatrix);
  double screenArray[4][4] = {
      {((double)width) / 2, 0, 0, ((double)width) / 2},
      {0, ((double)height) / 2, 0, ((double)height) / 2},
      {0, 0, 1, 0},
      {0, 0, 0, 1}};
  cameraMatrix = cameraMatrix.transpose();
  projMatrix = projMatrix.transpose(); // Cause maya representation
  screenMatrix = MMatrix(screenArray);

  MPoint controlScreenPt = screenMatrix * projMatrix * cameraMatrix * controlPt;
  controlScreenPt.cartesianize();
  offsetVec.x = controlScreenPt.x - curr_x;
  offsetVec.y = controlScreenPt.y - curr_y;
  offsetVec.z = controlScreenPt.z - 0;

  editDragMode = true;
  strokeColor = MColor(0.2, 1.0, 0.4, 0);
}

MStatus ControlInteraction::doRelease(MEvent &event,
                                      MHWRender::MUIDrawManager &drawMgr,
                                      const MHWRender::MFrameContext &context) {
  removeCallback(); // cout << "ControlInteraction::doRelease" << endl;

  event.getPosition(curr_x, curr_y);
  screenStrokePoints.append(MPoint(curr_x, curr_y, 0));
  last_x = curr_x;
  last_y = curr_y;

  if (screenStrokePoints.length() >=
      7) { // If not tapping and have enough points for a stroke.
    MGlobal::setActiveSelectionList(prevSelection);

    MPoint pathPoint = MPoint::origin;
    if (editDragMode) { // Dragging finished.  No need to set again.
      editDragMode = false;
      setUI();
    } else { // Handle Query mode or QM mode
      int handlePointsCode =
          squidgetToolManager->handlePoints(screenStrokePoints, pathPoint);

      if (handlePointsCode == 2) { // Stroke off canvas, use BookmarkSquidgetTool
        // if (queryMode) { // Query mode
        //   cout << "Query mode" << endl;
        //   tool->setPenControlPoints(screenStrokePoints);
        //   tool->queryClosestObjectOnly();
        // } else {
        //   tool->setPenControlPoints(screenStrokePoints);
        //   tool->predoIt();
        //   // if not exist otherwise use previously existing query
        //   // QueryClosestObjectOnly()
        //   tool->redoIt();
        // }
      }
    }
  }

  tool->finalize(); // Reset Everything to normal state.
  screenStrokePoints.clear();
  queryMode = false;
  squidgetToolManager->unHandleCanvas();
  return MS::kSuccess;
}

MStatus
ControlInteraction::doPtrMoved(MEvent &event,
                               MHWRender::MUIDrawManager &drawMgr,
                               const MHWRender::MFrameContext &context) {
  event.getPosition(curr_x, curr_y);
  return MS::kSuccess;
}

MStatus ControlInteraction::handleTool(BookmarkSquidgetTool *t) {
  tool = t;
  return MS::kSuccess;
}

bool ControlInteraction::completeAction() {
  return squidgetToolManager
      ->toggleSquidgetCanvasesVisibility(); // Toggle hidding canvases.
}

void ControlInteraction::deleteAction() {}

void ControlInteraction::abortAction() {}

void ControlInteraction::draw(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                              const MHWRender::MFrameContext &context) {
  draw(drawMgr);
}

void ControlInteraction::draw(MHWRender::MUIDrawManager &drawMgr) {
  drawMgr.beginDrawable();
  drawMgr.setColor(strokeColor);
  drawMgr.setPointSize((float)cursorSize);
  drawMgr.setLineWidth((float)strokeWidth);
  drawMgr.setLineStyle(MUIDrawManager::kSolid);
  drawMgr.lineStrip(screenStrokePoints, true);
  drawMgr.endDrawable();
  refresh();
}

INTERACTION_MODE ControlInteraction::getMode() {
  return INTERACTION_MODE::EDIT_MODE;
}

void ControlInteraction::setUI() {
  strokeColor = MColor(0, 0.7, 1);
  strokePointDist = 4;
  strokeWidth = 3;
  cursorSize = 3;
}

} // namespace interaction
} // namespace apps
} // namespace bookmark_squidget
