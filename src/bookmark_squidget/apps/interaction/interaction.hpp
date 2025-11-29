#ifndef INTERACTION_MODE_HPP
#define INTERACTION_MODE_HPP

#include <memory>

#include <maya/M3dView.h>
#include <maya/MEvent.h>
#include <maya/MHardwareRenderer.h>
#include <maya/MStatus.h>

#include "bookmark_squidget/canvas/squidgetToolManager.hpp"

namespace bookmark_squidget {
namespace apps {
namespace interaction {
enum INTERACTION_MODE {
  NONE,
  CREATE_MODE,
  EDIT_MODE,
  EDIT_DRAG_MODE,
  GLOBAL_MODE,
};

class InteractionMode {
public:
  virtual ~InteractionMode() {}

  virtual void modeSetup() = 0;
  virtual void modeCleanup() = 0;

  virtual MStatus doPress(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                          const MHWRender::MFrameContext &context) = 0;
  virtual MStatus doDrag(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                         const MHWRender::MFrameContext &context) = 0;
  virtual MStatus doRelease(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                            const MHWRender::MFrameContext &context) = 0;
  virtual MStatus doPtrMoved(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                             const MHWRender::MFrameContext &context) = 0;
  // virtual MStatus handleTool(std::shared_ptr<BookmarkSquidgetTool> tool) = 0;
  virtual MStatus handleTool(BookmarkSquidgetTool *tool) = 0;

  virtual bool completeAction() = 0;
  virtual void deleteAction() = 0;
  virtual void abortAction() = 0;

  virtual void draw(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                    const MHWRender::MFrameContext &context) = 0;
  virtual INTERACTION_MODE getMode() = 0;

  void refresh() {
    M3dView view = M3dView::active3dView();
    view.refresh();
  }

  std::shared_ptr<canvas::SquidgetToolManager> squidgetToolManager;
};
} // namespace interaction
} // namespace apps
} // namespace bookmark_squidget

#endif // INTERACTION_MODE_HPP