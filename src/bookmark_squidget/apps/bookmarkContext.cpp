#include "bookmark_squidget/apps/bookmarkContext.hpp"

namespace bookmark_squidget {
namespace apps {
BookmarkContext::BookmarkContext() {
  setTitleString("Bookmark Curve Tool");

  logger = std::make_shared<common::logging::Logger>("", "AbstractContext");
  curveHUD = std::make_shared<apps::CurveHUD>();
  squidgetToolManager = std::make_shared<canvas::SquidgetToolManager>();

  setCreateMode();
  // setEditMode();
  MGlobal::executePythonCommand("cmds.select(clear=True)");
}

BookmarkContext::~BookmarkContext() {}

void BookmarkContext::toolOnSetup(MEvent &event) {
  setHelpString("Draw a curve on the screen! ");

  setCreateMode();
  curveHUD->drawHUD();
  logger->logPluginActivated();
}

void BookmarkContext::toolOffCleanup() {
  interactionMode->modeCleanup();
  curveHUD->removeHUD();
  MPxContext::toolOffCleanup();
  logger->logPluginDeactivated();
}

// Viewport 2.0 Methods ########################################################
MStatus BookmarkContext::doPress(MEvent &event,
                                 MHWRender::MUIDrawManager &drawMgr,
                                 const MHWRender::MFrameContext &context) {
  mouseButton = event.mouseButton();
  event.getPosition(curr_x, curr_y);
  int mouseNum = (mouseButton & MEvent::kMiddleMouse) ? 1 : 0;
  logger->logMousePress(mouseNum, curr_x, curr_y);

  interactionMode->doPress(event, drawMgr, context);
  interactionMode->draw(event, drawMgr, context);
  if (interactionMode->getMode() == interaction::EDIT_MODE) {
    BookmarkSquidgetTool *tool = (BookmarkSquidgetTool *)newToolCommand();
    interactionMode->handleTool(tool);
  }
  return MS::kSuccess;
}

MStatus BookmarkContext::doDrag(MEvent &event,
                                MHWRender::MUIDrawManager &drawMgr,
                                const MHWRender::MFrameContext &context) {
  interactionMode->doDrag(event, drawMgr, context);
  interactionMode->draw(event, drawMgr, context);
  return MS::kSuccess;
}

MStatus BookmarkContext::doRelease(MEvent &event,
                                   MHWRender::MUIDrawManager &drawMgr,
                                   const MHWRender::MFrameContext &context) {
  mouseButton = event.mouseButton();
  modifier = event.modifiers();
  event.getPosition(curr_x, curr_y);

  int mouseNum = (mouseButton & MEvent::kMiddleMouse) ? 1 : 0;
  logger->logMouseRelease(mouseNum, curr_x, curr_y);

  MStatus status = interactionMode->doRelease(event, drawMgr, context);
  interactionMode->draw(event, drawMgr, context);
  return MS::kSuccess;
}

MStatus BookmarkContext::doPtrMoved(MEvent &event,
                                    MHWRender::MUIDrawManager &drawMgr,
                                    const MHWRender::MFrameContext &context) {
  interactionMode->doPtrMoved(event, drawMgr, context);
  interactionMode->draw(event, drawMgr, context);
  return MS::kSuccess;
}

MStatus BookmarkContext::drawFeedback(MHWRender::MUIDrawManager &drawMgr,
                                      const MHWRender::MFrameContext &context) {
  return MS::kSuccess;
}

// Keyboard Actions ### ########################################################
void BookmarkContext::completeAction() {
  logger->logCompleteAction();
  bool good = interactionMode->completeAction();

  M3dView view = M3dView::active3dView();
  view.refresh();
  curveHUD->setCanvasVisible(good);
  curveHUD->updateHUD();
}

void BookmarkContext::deleteAction() {
  logger->logDeleteAction();
  interactionMode->deleteAction();

  M3dView view = M3dView::active3dView();
  view.refresh();
}

void BookmarkContext::abortAction() {
  logger->logAbortAction();

  interactionMode->abortAction();
  toggleMode();

  M3dView view = M3dView::active3dView();
  view.refresh();
}

void BookmarkContext::toggleMode() {
  if (interactionMode->getMode() == interaction::CREATE_MODE) {
    setEditMode();
  } else if (interactionMode->getMode() == interaction::EDIT_MODE) {
    setCreateMode();
  }
  curveHUD->updateHUD();
  M3dView view = M3dView::active3dView();
  view.refresh();
}

void BookmarkContext::setCreateMode() {
  logger->logCreateMode();
  interactionMode =
      std::make_shared<interaction::CreateInteraction>(squidgetToolManager);
  MGlobal::displayInfo("Create Mode");
  curveHUD->setMode(true);
  setCursor(MCursor::defaultCursor);
}

void BookmarkContext::setEditMode() {
  logger->logEditMode();
  interactionMode =
      std::make_shared<interaction::ControlInteraction>(squidgetToolManager);
  MGlobal::displayInfo("Control Mode");
  curveHUD->setMode(false);
  setCursor(MCursor::pencilCursor);
}

} // namespace apps
} // namespace bookmark_squidget
