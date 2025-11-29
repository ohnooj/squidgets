#include "abstract_squidget/apps/abstractContext.hpp"

#include <maya/M3dView.h>
#include <maya/MFnTransform.h>
#include <maya/MGlobal.h>
#include <maya/MMatrix.h>
#include <maya/MTimerMessage.h>
#include <maya/MUIDrawManager.h>

#include "abstract_squidget/apps/abstractSquidgetTool.hpp"
#include "abstract_squidget/deformer/abstractDeformer.hpp"
#include "common/util/math.hpp"
#include "common/util/maya.hpp"

#include "common/maya/factory.hpp"

namespace abstract_squidget {
namespace apps {
AbstractContext::AbstractContext() {
  cout << "AbstractContext()" << endl;
  setTitleString("Abstract Curve Tool");
  MGlobal::executePythonCommand("cmds.select(clear=True)");

  logger = std::make_shared<common::logging::Logger>("", "AbstractContext");

  squidgetType = std::make_shared<MString>("Abstract");

  objCandidatesStrArray = MStringArray();
  objCandidatesStrArray.append("lemon2");

  editedObj = std::make_shared<MString>("");

  shadowCurveLocalStr = std::make_shared<MString>("");
  shadowCurveWorld = std::make_shared<MPointArray>(0);
  resetContextState(false);
}

AbstractContext::~AbstractContext() { cout << "~AbstractContext()" << endl; }

void AbstractContext::toolOnSetup(MEvent &event) {
  logger->logPluginActivated();
  setHelpString("Draw a curve on the screen! ");
  // MGlobal::executeCommand("source abstractUI.mel; abstractToolSettingsUI;");

  resetContextState(false);
}

void AbstractContext::toolOffCleanup() {
  MGlobal::executeCommand("if (`control -exists abstractToolSettingsForm`) "
                          "deleteUI abstractToolSettingsForm;");
  MPxContext::toolOffCleanup();
  logger->logPluginDeactivated();
}

void AbstractContext::getClassName(MString &name) const {
  name.set("abstractTool");
}

// =============================================================================
//                          Viewport 2.0 Methods
// =============================================================================
MStatus AbstractContext::doPress(MEvent &event,
                                 MHWRender::MUIDrawManager &drawMgr,
                                 const MHWRender::MFrameContext &context) {
  removeCallback();
  mouseButton = event.mouseButton();
  modifier = event.modifiers();
  preSelectMode = event.isModifierShift();
  deformMode = taskMode == 3;
  translateOnlyMode = taskMode == 1;

  if (preSelectMode) {
    strokeColor = SHIFT_COLOR;
    logger->logShift();
  } else if (deformMode) {
    strokeColor = CTRL_COLOR;
  } else if (translateOnlyMode) {
    strokeColor = MMB_COLOR;
  }

  if (!verifyInputMode() || !verifyTaskMode()) {
    draw(event, drawMgr, context);
    return MS::kFailure;
  }

  // Set initial curve point
  event.getPosition(curr_x, curr_y);
  screenStrokePoints.clear();
  addStrokePoint();

  // Log
  int mouseNum = (mouseButton & MEvent::kMiddleMouse) ? 1 : 0;
  logger->logMousePress(mouseNum, curr_x, curr_y);

  // common::util::log::log("doPress");
  MEvent eventCopy = event;
  if (taskMode == 2) {
    startCallback(curr_x, curr_y, event);
  }

  draw(event, drawMgr, context);
  return MS::kSuccess;
}

MStatus AbstractContext::doDrag(MEvent &event,
                                MHWRender::MUIDrawManager &drawMgr,
                                const MHWRender::MFrameContext &context) {
  mouseButton = event.mouseButton();
  modifier = event.modifiers();
  if (!verifyInputMode() || !verifyTaskMode()) {
    draw(event, drawMgr, context);
    return MS::kFailure;
  }
  // cout << "AbstractContext::doDrag()" << endl;

  short x, y;
  event.getPosition(x, y);
  MStatus status = doDrag(x, y, event);

  draw(event, drawMgr, context);
  return MS::kSuccess;
}

MStatus AbstractContext::doDrag(short x, short y, MEvent &event) {
  curr_x = x;
  curr_y = y;

  double strokeFarEnough =
      sqrt(pow(curr_x - last_x, 2) + pow(curr_y - last_y, 2)) > strokePointDist;
  if (strokeFarEnough) {
    addStrokePoint();
  }

  if (!editDragMode) {
    startCallback(curr_x, curr_y, event);
  }

  if (editDragMode) {
    removeCallback();
    shadowCurveWorld->clear();
    // cout << "AbstractContext::doDrag() editDragMode" << endl;
    MPoint mouse(curr_x, curr_y, 0);

    bool mmb = event.mouseButton() & MEvent::kMiddleMouse;
    bool canMod = (!translateOnlyMode) && event.isModifierControl();

    offsetVec = deformer::livePoint(*editedObj, mouse, offsetVec, canMod);
  }
  return MS::kSuccess;
}

MStatus AbstractContext::doRelease(MEvent &event,
                                   MHWRender::MUIDrawManager &drawMgr,
                                   const MHWRender::MFrameContext &context) {
  MStatus status;
  mouseButton = event.mouseButton();
  modifier = event.modifiers();

  event.getPosition(curr_x, curr_y);
  addStrokePoint();
  removeCallback();

  // Log
  int mouseNum = (mouseButton & MEvent::kMiddleMouse) ? 1 : 0;
  logger->logMouseRelease(mouseNum, curr_x, curr_y);

  // If is tap, then select the object
  bool longEnoughStroke = screenStrokePoints.length() > 4;
  if (!longEnoughStroke && !editDragMode) {
    // select the object if the stroke is too short.

    logger->logMouseTap();

    status = MGlobal::selectFromScreen(screenStrokePoints[0].x,
                                       screenStrokePoints[0].y,
                                       MGlobal::kReplaceList);
    resetContextState(false);
    draw(event, drawMgr, context);
    return MS::kSuccess;
  }

  if (!verifyInputMode() || !verifyTaskMode()) {
    draw(event, drawMgr, context);
    return MS::kFailure;
  }

  if (editDragMode) {
    // Finish drag mode
    screenStrokePoints.clear(); // Pass an empty array to do nothing
  } else {
    bool mmb = event.mouseButton(&status) == MEvent::kMiddleMouse;
    bool isRigCurves = *squidgetType == "Rig";
    CallToolArguments callToolContextArgs(
        getStrokePointsStr(), getObjCandidatesStr(), preSelectMode, deformMode,
        translateOnlyMode, *shadowCurveLocalStr, isRigCurves);

    callToolContext(callToolContextArgs, status);
  }

  resetContextState(preSelectMode || isInputMode(3));
  draw(event, drawMgr, context);
  return status;
}

MStatus AbstractContext::doPtrMoved(MEvent &event,
                                    MHWRender::MUIDrawManager &drawMgr,
                                    const MHWRender::MFrameContext &context) {
  // cout << "AbstractContext::doPtrMoved()" << endl;
  event.getPosition(curr_x, curr_y);
  draw(event, drawMgr, context);
  return MS::kSuccess;
}

MStatus AbstractContext::drawFeedback(MHWRender::MUIDrawManager &drawMgr,
                                      const MHWRender::MFrameContext &context) {
  draw(drawMgr);
  return MS::kSuccess;
}

// =============================================================================
//                            Keyboard Actions
// =============================================================================
void AbstractContext::completeAction() {
  logger->logCompleteAction();
  M3dView view = M3dView::active3dView();
  view.refresh();
}

void AbstractContext::deleteAction() {
  logger->logDeleteAction();
  M3dView view = M3dView::active3dView();
  view.refresh();
}

void AbstractContext::abortAction() {
  logger->logAbortAction();
  resetContextState(false);
  M3dView view = M3dView::active3dView();
  view.refresh();
}

// =============================================================================
//                              CALLBACKS
// =============================================================================
void AbstractContext::startCallback(short x, short y, MEvent &event) {
  if (inputMode == 1 || inputMode == 3) {
    return;
  }

  // The HOLD functionality for MPxContext doesn't work as expected.  We need
  // need to use a callback timer when the mouse is held down and stationary.
  // We have one callback at a time, and if the mouse moves, then we remove the
  // callback and set a new one.
  // If the callback lives after .4 seconds, then we perform a regular
  // 1 stroke AbstractSquidgetTool action to the curve, but then we transform
  // the object around OUTSIDE THE TOOL.  So there is not maya command being
  // called when I move the live.
  // When we release our mouse, we call a 1 stroke AbstractSquidgetTool but with
  // stripped down arguments.
  removeCallback();
  MyTimerCallbackData *callbackData =
      new MyTimerCallbackData(this, x, y, std::make_shared<MEvent>(event));
  callbackId =
      MTimerMessage::addTimerCallback(.4f, timerCallback, callbackData);
}

void AbstractContext::removeCallback() {
  // If a callback is already registered, remove it.
  if (callbackId != static_cast<MCallbackId>(0)) {
    MTimerMessage::removeCallback(callbackId);
  }
  callbackId = static_cast<MCallbackId>(0);
}

void AbstractContext::timerCallback(float elapsedTime, float lastTime,
                                    void *clientData) {
  MyTimerCallbackData *obj = static_cast<MyTimerCallbackData *>(clientData);
  if (obj && obj->obj && obj->event && !obj->obj->isDragMode()) {
    if (obj->obj->isDragMode()) {
      // cout << "Why am i here" << endl;
      return;
    }
    if (obj->obj->screenStrokePoints.length() < 5) {
      // cout << "Too short of a stroke" << endl;
      return;
    }
    obj->obj->removeCallback();

    obj->obj->logger->logMouseHold();

    obj->obj->shadowCurveWorld->clear();

    obj->obj->setDragMode();
    obj->obj->doDrag(obj->x, obj->y, *(obj->event));

    M3dView view = M3dView::active3dView();
    view.refresh();
  } else {
    // cout << "timerCallback() failed" << endl;
  }
}

bool AbstractContext::isDragMode() { return editDragMode; }

void AbstractContext::setDragMode() {
  MPoint controlPt = MPoint::origin;

  // Call the tool context to
  MStatus status;

  // bool mmb = event.mouseButton(&status) & MEvent::kMiddleMouse;
  bool isRigCurves = *squidgetType == "Rig";
  CallToolArguments callToolContextArgs(
      getStrokePointsStr(), getObjCandidatesStr(), preSelectMode, deformMode,
      translateOnlyMode, "", isRigCurves);
  callToolContext(callToolContextArgs, status);
  // shadowCurvesStr = ""; // remove this for later

  MDagPath geomPath = common::util::getDagPathFromName(*editedObj);
  MFnTransform fnTransform(geomPath);
  controlPt = fnTransform.getTranslation(MSpace::kWorld, &status);

  MPoint controlScreenPt = common::util::getWorldToScreenPoint(controlPt);
  controlScreenPt.cartesianize();
  offsetVec.x = curr_x - controlScreenPt.x;
  offsetVec.y = curr_y - controlScreenPt.y;
  offsetVec.z = 0 - controlScreenPt.z;

  editDragMode = true;
  strokeColor = transparent_color; // Hide Stroke
  screenStrokePoints.clear();
  view = M3dView::active3dView();
  view.refresh();

  editDragMode = true;
}

// =============================================================================
//                              TOOL ACTIONS
// =============================================================================
MString AbstractContext::callToolContext(CallToolArguments args,
                                         MStatus &status) {
  AbstractSquidgetTool *tool = (AbstractSquidgetTool *)newToolCommand();

  // Setup Arguments
  tool->setCurvePoints2D(args.strokePointsStr); // Pass in empty to do nothing
  tool->setCandidateObjs(args.objCandidatesStr);

  if (args.deform) {
    tool->setDeform(args.deform);
  } else if (args.translateOnly) {
    tool->setTranslateOnly(args.translateOnly);
  }

  if (args.preSelect) {
    tool->setPreSelect(args.preSelect); // Preselect and object
  } else if (args.curveToMatchStr != "") {
    // If we are in the 2 stroke mode, then we need to pass in preselect
    tool->setCurveToMatchStr(args.curveToMatchStr);
  }

  tool->setIsRigCurves(args.isRigCurves);

  tool->fromContextShadowCurveWorld(shadowCurveWorld);
  tool->fromContextShadowCurveLocalStr(shadowCurveLocalStr);
  tool->fromContextEditedObject(editedObj);
  // tool->set

  tool->setLogger(logger);
  logger->logSketchStroke(args.strokePointsStr.asChar());

  // Call tool
  // Automatically updates the contexts
  tool->redoIt();

  logger->logCommand(tool->getCallableCommandStr().asChar());
  logger->logSelectedObjectName(editedObj->asChar());

  // ConvertshadowCurveWorld to string for LOGGING
  MString shadowCurveWorldStr = "";
  for (unsigned int i = 0; i < shadowCurveWorld->length(); i++) {
    MPoint shadowPt = (*shadowCurveWorld)[i];
    short x, y;
    M3dView view = M3dView::active3dView();
    view.worldToView(shadowPt, x, y);

    shadowCurveWorldStr += x;
    shadowCurveWorldStr += " ";
    shadowCurveWorldStr += y;
    shadowCurveWorldStr += " ";
    shadowCurveWorldStr += shadowPt.z;
    if (i != shadowCurveWorldStr.length() - 1) {
      shadowCurveWorldStr += " ";
    }
  }
  logger->logProjectedCurve(shadowCurveWorldStr.asChar());

  // Finalize call
  status = tool->finalize();

  MGlobal::selectByName(*editedObj, MGlobal::kReplaceList);
  return *editedObj;
}

// =============================================================================
//                          HELPER METHODS
// =============================================================================
bool AbstractContext::verifyInputMode() {
  /**
   * @brief
   * 1 Stroke: draw a stroke and no shift
   * 1 Stroke Live: draw a stroke and hold mouse (ctrl only if editDragMode)
   * 2 Stroke: draw a stroke + <shift> -> draw another stroke.
   */
  bool shift = modifier & MEvent::shiftKey;

  bool free_form = inputMode == 0;
  bool one_stroke = inputMode == 1 && !shift;
  bool one_stroke_live = inputMode == 2 && !shift;

  bool isFirstStroke = shift;
  bool isSecondStroke = !shift && *shadowCurveLocalStr != "";
  bool two_stroke = inputMode == 3 && (isFirstStroke || isSecondStroke);

  return free_form || one_stroke || one_stroke_live || two_stroke;
}

bool AbstractContext::isInputMode(int mode) { return inputMode == mode; }

bool AbstractContext::verifyTaskMode() {
  /**
   * @brief
   * Translate: Middle Mouse
   * Translate + Rotate: Left Mouse
   * Deform:
   */
  bool left_mouse = mouseButton & MEvent::kLeftMouse;
  bool middle_mouse = mouseButton & MEvent::kMiddleMouse;
  bool ctrl = modifier & MEvent::controlKey;

  bool freeTask = taskMode == 0;
  bool translateTask = taskMode == 1 && (!middle_mouse && left_mouse) && !ctrl;
  bool translateRotTask = taskMode == 2 && (!middle_mouse && left_mouse) &&
                          ((!editDragMode && !ctrl) || (editDragMode));
  bool deformTask = taskMode == 3 && (!middle_mouse && left_mouse);

  return freeTask || translateTask || translateRotTask || deformTask;
}

bool AbstractContext::isTaskMode(int mode) { return taskMode == mode; }

bool AbstractContext::resetContextState(bool keepPreviousToolResults) {
  screenStrokePoints.clear();
  callbackId = -1;
  setUI();

  editDragMode = false;
  preSelectMode = false;
  deformMode = false;
  if (!keepPreviousToolResults) {
    editedObj->clear();
    shadowCurveLocalStr->clear();
    shadowCurveWorld->clear();
  }

  return true;
}

bool AbstractContext::addStrokePoint() {
  return addStrokePoint(curr_x, curr_y);
}

bool AbstractContext::addStrokePoint(short x, short y) {
  MPoint newPoint = MPoint(x, y, 0);
  screenStrokePoints.append(newPoint);
  last_x = x;
  last_y = y;
  return true;
}

// =============================================================================
//                              Draw Methods
// =============================================================================
void AbstractContext::draw(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                           const MHWRender::MFrameContext &context) {
  draw(drawMgr);
}

void AbstractContext::draw(MHWRender::MUIDrawManager &drawMgr) {
  drawMgr.beginDrawable();
  drawMgr.setColor(strokeColor);
  drawMgr.setPointSize((float)cursorSize);
  drawMgr.setLineWidth((float)strokeWidth);
  drawMgr.setLineStyle(MUIDrawManager::kSolid);
  drawMgr.lineStrip(screenStrokePoints, true);
  drawMgr.endDrawable();

  displayShadowCurveWorld(drawMgr);
  view = M3dView::active3dView();
  view.refresh();
}

void AbstractContext::displayShadowCurveWorld(
    MHWRender::MUIDrawManager &drawMgr) {
  if (shadowCurveWorld->length() == 0 || inputMode != 3) {
    return;
  }

  MPointArray shadowCurveScreen;
  M3dView view = M3dView::active3dView();
  for (unsigned int i = 0; i < shadowCurveWorld->length(); i++) {
    MPoint shadowPt = (*shadowCurveWorld)[i];
    MPoint shadowScreenPt = common::util::getWorldToScreenPoint(shadowPt);
    shadowCurveScreen.append(shadowScreenPt);
  }

  drawMgr.beginDrawable();
  drawMgr.setColor(shadowColor);
  drawMgr.setPointSize((float)cursorSize);
  drawMgr.setLineWidth((float)strokeWidth);
  drawMgr.setLineStyle(MUIDrawManager::kSolid);
  drawMgr.lineStrip(shadowCurveScreen, true);
  drawMgr.endDrawable();
}

void AbstractContext::setUI() { strokeColor = LMB_COLOR; }

// =============================================================================
// Getters and Setters
// =============================================================================
void AbstractContext::setStrokeWidth(double width) {
  strokeWidth = width;
  resetContextState(false);
}

int AbstractContext::getStrokeWidth() { return strokeWidth; }

void AbstractContext::setStrokePointDist(double dist) {
  strokePointDist = dist;
  resetContextState(false);
}

int AbstractContext::getStrokePointDist() { return strokePointDist; }

void AbstractContext::setInputMode(int mode) {
  this->inputMode = mode;
  logger->logInputMode(inputMode);
  resetContextState(false);
}

int AbstractContext::getInputMode() { return inputMode; }

void AbstractContext::setTaskMode(int mode) {
  this->taskMode = mode;
  logger->logTaskMode(taskMode);
  resetContextState(false);
}

int AbstractContext::getTaskMode() { return taskMode; }

void AbstractContext::setObjCandidatesStr(MString objCandidatesStr) {
  // Split objCandidates by space and add to objCandidates array.
  objCandidatesStrArray.clear();
  objCandidatesStr.split(' ', objCandidatesStrArray);
  resetContextState(false);
}

MString AbstractContext::getObjCandidatesStr() {
  if (*editedObj != "" && !preSelectMode) {
    // I preselected and object and want to change my preseleciton.
    return *editedObj;
  }

  // Combine objCandidates array into a single string with spaces.
  MString ret = "";
  for (unsigned int i = 0; i < objCandidatesStrArray.length(); i++) {
    ret += objCandidatesStrArray[i];
    if (i != objCandidatesStrArray.length() - 1) {
      ret += " ";
    }
  }

  return ret;
}

void AbstractContext::setLoggerPath(MString loggerPath) {
  logger->setLogFilepath(loggerPath.asChar());
  resetContextState(false);
}

MString AbstractContext::getLoggerPath() {
  return MString(logger->getLogFilepath().c_str());
}

void AbstractContext::setSquidgetType(MString squidgetType) {
  MStatus status;
  if (squidgetType == "Abstract") {
    this->squidgetType = std::make_shared<MString>("Abstract");
  } else if (squidgetType == "Rig") {
    this->squidgetType = std::make_shared<MString>("Rig");
  } else {
    status = MStatus::kFailure;
    CHECK_MSTATUS(status);
  }
  logger->logSquidgetType(squidgetType.asChar());

  this->squidgetType = std::make_shared<MString>(squidgetType);
  resetContextState(false);
}

MString AbstractContext::getSquidgetType() { return *squidgetType; }

MString AbstractContext::getStrokePointsStr() {
  MString ret = "";
  ret += screenStrokePoints[0].x;
  ret += " ";
  ret += screenStrokePoints[0].y;
  ret += " ";
  ret += screenStrokePoints[0].z;
  for (unsigned int i = 1; i < screenStrokePoints.length(); i++) {
    ret += " ";
    ret += screenStrokePoints[i].x;
    ret += " ";
    ret += screenStrokePoints[i].y;
    ret += " ";
    ret += screenStrokePoints[i].z;
  }

  return ret;
}

} // namespace apps
} // namespace abstract_squidget
