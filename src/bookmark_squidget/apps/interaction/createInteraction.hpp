#ifndef CREATE_INTERACTION_HPP
#define CREATE_INTERACTION_HPP

#include <maya/MEvent.h>
#include <maya/MHardwareRenderer.h>
#include <maya/MStatus.h>
#include <maya/MUIDrawManager.h>

#include "bookmark_squidget/apps/interaction/interaction.hpp"
#include "bookmark_squidget/canvas/squidgetToolManager.hpp"

namespace bookmark_squidget {
namespace apps {
namespace interaction {

class CreateInteraction : public InteractionMode {
public:
  CreateInteraction(
      std::shared_ptr<canvas::SquidgetToolManager> squidgetToolManager);
  ~CreateInteraction() {}

  void modeSetup();
  void modeCleanup();

  MStatus doPress(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                  const MHWRender::MFrameContext &context);
  MStatus doDrag(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                 const MHWRender::MFrameContext &context);
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
  INTERACTION_MODE getMode();

private:
  void setUI();

  MStatus handleSelections(MEvent &event);

  MColor strokeColor;
  double strokePointDist;
  double strokeWidth;
  double cursorSize;

  MPointArray screenStrokePoints;
  MSelectionList prevSelection;
  short curr_x, curr_y;
  short last_x, last_y;
};
} // namespace interaction
} // namespace apps
} // namespace bookmark_squidget

#endif // CREATE_INTERACTION_HPP