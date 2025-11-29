#ifndef EDIT_INTERACTION_HPP
#define EDIT_INTERACTION_HPP

#include <maya/MDGMessage.h>
#include <maya/MEvent.h>
#include <maya/MHardwareRenderer.h>
#include <maya/MStatus.h>
#include <maya/MTimerMessage.h>
#include <maya/MUIDrawManager.h>

// This needs to come before the next include for some reason?
#include "common/util/math.hpp"

#include "bookmark_squidget/apps/interaction/interaction.hpp"
#include "bookmark_squidget/canvas/squidgetToolManager.hpp"

namespace bookmark_squidget {
namespace apps {
namespace interaction {

// struct callbackData {
//   bool* editDragMode;
//   MColor* strokeColor;
// };

class ControlInteraction : public InteractionMode {
public:
  ControlInteraction(
      std::shared_ptr<canvas::SquidgetToolManager> squidgetToolManager);
  ~ControlInteraction() {}

  void modeSetup();
  void modeCleanup();

  MStatus doPress(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                  const MHWRender::MFrameContext &context);
  MStatus doDrag(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                 const MHWRender::MFrameContext &context);
  // MStatus doDrag(MEvent &event, MHWRender::MUIDrawManager &drawMgr);
  MStatus doDrag(short x, short y, bool isCtrl,
                 MHWRender::MUIDrawManager &drawMgr);
  MStatus doRelease(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                    const MHWRender::MFrameContext &context);
  MStatus doPtrMoved(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                     const MHWRender::MFrameContext &context);
  // MStatus handleTool(std::shared_ptr<BookmarkSquidgetTool> tool);
  MStatus handleTool(BookmarkSquidgetTool *tool);

  bool completeAction();
  void deleteAction();
  void abortAction();

  void draw(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
            const MHWRender::MFrameContext &context);
  void draw(MHWRender::MUIDrawManager &drawMgr);

  INTERACTION_MODE getMode();
  void setDragMode();
  bool isDragMode();

private:
  void setUI();

  MCallbackId callbackId = -1;
  void startCallback(short x, short y, MHWRender::MUIDrawManager *drawMgr);
  void removeCallback();
  static void timerCallback(float elapsedTime, float lastTime,
                            void *clientData);

  MColor strokeColor;
  double strokePointDist;
  double strokeWidth;
  double cursorSize;

  MPointArray screenStrokePoints;
  MSelectionList prevSelection;
  short curr_x, curr_y;
  short last_x, last_y;

  std::chrono::milliseconds dragTime;
  bool editDragMode;
  bool queryMode;
  MVector offsetVec;

  BookmarkSquidgetTool *tool;

  struct MyTimerCallbackData {
    ControlInteraction *obj;
    short x, y;
    MUIDrawManager *uiDrawManager;
    // MyTimerCallbackData(ControlInteraction* o, MEvent* e, MUIDrawManager*
    // udm)
    //         : obj(o), event(e), uiDrawManager(udm){}
    MyTimerCallbackData(ControlInteraction *o, short x, short y,
                        MUIDrawManager *udm)
        : obj(o), x(x), y(y), uiDrawManager(udm) {}
  };
};

} // namespace interaction
} // namespace apps
} // namespace bookmark_squidget

#endif // EDIT_INTERACTION_HPP
