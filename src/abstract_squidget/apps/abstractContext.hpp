#ifndef ABSTRACT_CONTEXT_HPP
#define ABSTRACT_CONTEXT_HPP

#include <memory>

#include <maya/M3dView.h>
#include <maya/MCallbackIdArray.h>
#include <maya/MColor.h>
#include <maya/MPointArray.h>
#include <maya/MPxSelectionContext.h>
#include <maya/MSelectionList.h>
#include <maya/MStringArray.h>

#include "common/logging/logger.hpp"

/**
 * @brief Context that draws 2D strokes on the viewport 2.0.  This class mostly
 * relies on MHWRender::MUIDrawManager within the mouse methods to display the
 * curve.  There are simpler mouse input methods that rely on OpenGL, but I do
 * not implement these.
 *
 * This class should:
 *  - draw graphics on the viewport using I/O input.
 *  - pass mouse inputs into CurveTool.
 *
 * This class should not:
 *  - Deform geometry.  Logic should be handled in CurveTool.
 *
 * Handles everything in 2D screen space.
 */
namespace abstract_squidget {
namespace apps {
class AbstractContext : public MPxContext {

  struct MyTimerCallbackData {
    AbstractContext *obj;
    short x, y;
    std::shared_ptr<MEvent> event;
    // MHWRender::MUIDrawManager *uiDrawManager;
    // MyTimerCallbackData(ControlInteraction* o, MEvent* e, MUIDrawManager*
    // udm)
    //         : obj(o), event(e), uiDrawManager(udm){}
    MyTimerCallbackData(AbstractContext *o, short x, short y,
                        std::shared_ptr<MEvent> event)
        // MHWRender::MUIDrawManager *udm)
        : obj(o), x(x), y(y), event(event) {}
    // : obj(o), x(x), y(y), event(event), uiDrawManager(udm) {}
  };

  struct CallToolArguments {
    MString strokePointsStr;
    MString objCandidatesStr;
    bool preSelect = false;
    bool deform = false;
    bool translateOnly = false;
    MString curveToMatchStr;
    bool isRigCurves = false;

    CallToolArguments(MString strokePointsStr, MString objCandidatesStr,
                      bool preSelect, bool deform, bool translateOnly,
                      MString shadowCurveStr, bool isRigCurves)
        : strokePointsStr(strokePointsStr), objCandidatesStr(objCandidatesStr),
          preSelect(preSelect), deform(deform), translateOnly(translateOnly),
          curveToMatchStr(shadowCurveStr), isRigCurves(isRigCurves) {}
  };

public:
  AbstractContext();
  ~AbstractContext();
  void toolOnSetup(MEvent &event) override;
  void toolOffCleanup() override;
  void getClassName(MString &name) const override;

  // Viewport 2.0 methods, will only be triggered in viewport 2.0.
  MStatus doPress(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                  const MHWRender::MFrameContext &context) override;
  MStatus doRelease(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                    const MHWRender::MFrameContext &context) override;
  MStatus doDrag(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                 const MHWRender::MFrameContext &context) override;
  MStatus doPtrMoved(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                     const MHWRender::MFrameContext &context) override;
  MStatus drawFeedback(MHWRender::MUIDrawManager &drawMgr,
                       const MHWRender::MFrameContext &context) override;

  void draw(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
            const MHWRender::MFrameContext &context);

  void completeAction() override;
  void deleteAction() override;
  void abortAction() override;

  void setStrokeWidth(double width);
  int getStrokeWidth();
  void setStrokePointDist(double dist);
  int getStrokePointDist();
  void setInputMode(int mode);
  int getInputMode();
  void setTaskMode(int mode);
  int getTaskMode();
  void setObjCandidatesStr(MString objCandidatesStr);
  MString getObjCandidatesStr();
  void setLoggerPath(MString loggerPath);
  MString getLoggerPath();
  void setSquidgetType(MString squidgetType);
  MString getSquidgetType();

  MString getStrokePointsStr();

  void setDragMode();
  bool isDragMode();

private:
  MStatus doDrag(short x, short y, MEvent &event);
  //  MHWRender::MUIDrawManager &drawMgr);
  void startCallback(short x, short y, MEvent &event);
  //  MHWRender::MUIDrawManager *drawMgr);
  void removeCallback();
  static void timerCallback(float elapsedTime, float lastTime,
                            void *clientData);
  void draw(MHWRender::MUIDrawManager &drawMgr);
  bool addStrokePoint();
  bool addStrokePoint(short x, short y);
  bool resetContextState(bool keepPreviousToolResults);
  MString callToolContext(CallToolArguments args, MStatus &status);

  // bool isModeInteraction(int mode);
  void displayShadowCurveWorld(MHWRender::MUIDrawManager &drawMgr);

  MColor strokeColor;
  double strokePointDist = 4;
  double strokeWidth = 3;
  double cursorSize = 3;

  MPointArray screenStrokePoints;
  MSelectionList prevSelection;
  short curr_x, curr_y;
  short last_x, last_y;

  MCallbackId callbackId = -1;
  bool editDragMode;
  bool preSelectMode;
  bool deformMode;
  bool translateOnlyMode;

  MEvent::MouseButtonType mouseButton;
  MEvent::ModifierType modifier;

  MStringArray objCandidatesStrArray; // List of object candidates
  // MString shadowCurvesStr; // Shadow curve in screen space
  std::shared_ptr<MString> editedObj;
  std::shared_ptr<MString> shadowCurveLocalStr;
  std::shared_ptr<MPointArray> shadowCurveWorld; // For displaying shadow curve

  MVector offsetVec;

  std::shared_ptr<common::logging::Logger> logger;

  bool verifyInputMode();
  bool verifyTaskMode();
  bool isInputMode(int mode);
  bool isTaskMode(int mode);
  int inputMode = 0;
  int taskMode = 0;

  // 'Abstract', 'Rig'
  // Not implemented: 'Bookmark'
  std::shared_ptr<MString> squidgetType;

  /*
  InputModes:
  0 - free
  1 - 1 Stroke
  2 - 1 Stroke + live
  3 - 2 Stroke

  TaskModes:
  0 - free
  1 - translate
  2 - T + R
  3 - Deform
  */

  MColor LMB_COLOR = MColor(0.129, 0.812, 1.0);
  MColor MMB_COLOR = MColor(0.98, 0.649, 0.118);
  MColor SHIFT_COLOR = MColor(0.406, 0.965, 0.5, 1.0);
  MColor transparent_color = MColor(0.0, 0.0, 0.0, 0.0);
  MColor shadowColor = MColor(0.96, 0.25, 0.25, 1.0);
  MColor CTRL_COLOR = MColor(0.863, 0.469, 0.918, 1.0);

  M3dView view;
  void setUI();
};
} // namespace apps
} // namespace abstract_squidget

#endif