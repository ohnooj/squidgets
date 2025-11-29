#include "bookmark_squidget/apps/interaction/createInteraction.hpp"

namespace bookmark_squidget {
namespace apps {
namespace interaction {

CreateInteraction::CreateInteraction(
    std::shared_ptr<canvas::SquidgetToolManager> squidgetToolManager) {
  // cout << "CreateInteraction()" << endl;
  this->squidgetToolManager = squidgetToolManager;
  modeSetup();
  // cout << "CreateInteraction() end" << endl;
}

void CreateInteraction::modeSetup() {
  setUI();
  squidgetToolManager->createMode();
}

void CreateInteraction::modeCleanup() {
  squidgetToolManager->hoverCanvas(-1, -1);
}

MStatus CreateInteraction::doPress(MEvent &event,
                                   MHWRender::MUIDrawManager &drawMgr,
                                   const MHWRender::MFrameContext &context) {
  MGlobal::getActiveSelectionList(prevSelection);
  event.getPosition(curr_x, curr_y);

  screenStrokePoints.clear();
  screenStrokePoints.append(MPoint(curr_x, curr_y, 0));

  last_x = curr_x;
  last_y = curr_y;
  return MS::kSuccess;
}

MStatus CreateInteraction::doDrag(MEvent &event,
                                  MHWRender::MUIDrawManager &drawMgr,
                                  const MHWRender::MFrameContext &context) {
  event.getPosition(curr_x, curr_y);

  // Append new point to curve if curr position is far enough from the last
  // point
  double distFromLast = sqrt(pow(curr_x - last_x, 2) + pow(curr_y - last_y, 2));
  if (distFromLast > strokePointDist) {
    screenStrokePoints.append(MPoint(curr_x, curr_y, 0));
    last_x = curr_x;
    last_y = curr_y;
  }
  return MS::kSuccess;
}

MStatus CreateInteraction::doRelease(MEvent &event,
                                     MHWRender::MUIDrawManager &drawMgr,
                                     const MHWRender::MFrameContext &context) {
  event.getPosition(curr_x, curr_y);
  screenStrokePoints.append(MPoint(curr_x, curr_y, 0));

  last_x = curr_x;
  last_y = curr_y;

  if (screenStrokePoints.length() < 4) { // Tapping
    screenStrokePoints.clear();
    return handleSelections(event);
  }

  MGlobal::setActiveSelectionList(prevSelection);
  MPoint pathPoint;
  int handlePointsCode =
      squidgetToolManager->handlePoints(screenStrokePoints, pathPoint);

  screenStrokePoints.clear();
  squidgetToolManager->unHandleCanvas();
  MGlobal::setActiveSelectionList(prevSelection);

  return MStatus::kSuccess;
}

MStatus CreateInteraction::doPtrMoved(MEvent &event,
                                      MHWRender::MUIDrawManager &drawMgr,
                                      const MHWRender::MFrameContext &context) {
  event.getPosition(curr_x, curr_y);
  squidgetToolManager->hoverCanvas(curr_x, curr_y);
  last_x = curr_x;
  last_y = curr_y;

  return MS::kSuccess;
}

MStatus CreateInteraction::handleTool(BookmarkSquidgetTool *tool) {
  // cout << "CreateInteraction::handleTool" << endl;
  tool->finalize();
  return MS::kSuccess;
}

bool CreateInteraction::completeAction() {
  squidgetToolManager->createSquidgetCanvas();
  return true;
}

void CreateInteraction::deleteAction() {}

void CreateInteraction::abortAction() {}

void CreateInteraction::draw(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                             const MHWRender::MFrameContext &context) {
  drawMgr.beginDrawable();
  drawMgr.setColor(strokeColor);
  drawMgr.setPointSize((float)cursorSize);
  drawMgr.setLineWidth((float)strokeWidth);
  drawMgr.setLineStyle(MUIDrawManager::kSolid);
  drawMgr.lineStrip(screenStrokePoints, true);
  drawMgr.endDrawable();

  refresh();
}

MStatus CreateInteraction::handleSelections(MEvent &event) {
  // MGlobal::getActiveSelectionList(prevSelection);
  short x, y;
  event.getPosition(x, y);

  if (event.modifiers() & MEvent::shiftKey) {
    MGlobal::selectFromScreen(x, y, MGlobal::kXORWithList);
  } else {
    MGlobal::selectFromScreen(x, y, MGlobal::kReplaceList);
  }

  return MS::kSuccess;
}

void CreateInteraction::setUI() {
  // cout << "CreateInteraction::setUI()" << endl;
  strokeColor = MColor(1, 0.5, 0.3);
  strokePointDist = 6;
  strokeWidth = 4;
  cursorSize = 5;
}

INTERACTION_MODE CreateInteraction::getMode() {
  return INTERACTION_MODE::CREATE_MODE;
}
} // namespace interaction
} // namespace apps
} // namespace bookmark_squidget