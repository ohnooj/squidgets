#ifndef BookmarkContext_HPP
#define BookmarkContext_HPP

#include <chrono>
#include <math.h>
#include <ostream>
#include <string>

#include <maya/MFnCamera.h>
#include <maya/MIOStream.h>

#include <maya/MEvent.h>
// #include <maya/MPxContext.h>
#include <maya/MPxSelectionContext.h>

#include <maya/M3dView.h>
#include <maya/MCursor.h>
#include <maya/MEvent.h>
#include <maya/MGlobal.h>
#include <maya/MPointArray.h>

#include <maya/MColor.h>
#include <maya/MFrameContext.h>
#include <maya/MPoint.h>
#include <maya/MToolsInfo.h>
#include <maya/MUIDrawManager.h>

#include <maya/MSpinLock.h>
#include <maya/MThreadAsync.h>

#include <maya/MFnNurbsSurface.h>
#include <maya/MGlobal.h>

#include "bookmark_squidget/apps/bookmarkSquidgetTool.hpp"
#include "bookmark_squidget/apps/curveHUD.hpp"
#include "bookmark_squidget/apps/interaction/controlInteraction.hpp"
#include "bookmark_squidget/apps/interaction/createInteraction.hpp"
#include "bookmark_squidget/apps/interaction/interaction.hpp"
#include "bookmark_squidget/canvas/squidget.hpp"
#include "bookmark_squidget/canvas/squidgetToolManager.hpp"
#include "common/logging/logger.hpp"

/**
 * @brief Context that draws 2D strokes on the viewport 2.0.  This class mostly
 * relies on MHWRender::MUIDrawManager within the mouse methods to display the
 * curve.  There are simpler mouse input methods that rely on OpenGL, but I do
 * not implement these.
 *
 * This class should:
 *  - draw graphics on the viewport using I/O input.
 *  - pass mouse inputs into BookmarkSquidgetTool.
 *
 * This class should not:
 *  - Deform geometry.  Logic should be handled in BookmarkSquidgetTool.
 *
 * Handles everything in 2D screen space.
 *
 */

namespace bookmark_squidget {
namespace apps {
class BookmarkContext : public MPxContext {
public:
  BookmarkContext();
  ~BookmarkContext();
  void toolOnSetup(MEvent &event) override;
  void toolOffCleanup() override;

  // Viewport 2.0 methods, will only be triggered in viewport 2.0.
  MStatus doPress(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                  const MHWRender::MFrameContext &context) override;
  MStatus doDrag(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                 const MHWRender::MFrameContext &context) override;
  MStatus doRelease(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                    const MHWRender::MFrameContext &context) override;
  MStatus doPtrMoved(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                     const MHWRender::MFrameContext &context) override;
  MStatus drawFeedback(MHWRender::MUIDrawManager &drawMgr,
                       const MHWRender::MFrameContext &context) override;

  void completeAction() override;
  void deleteAction() override;
  void abortAction() override;

private:
  void toggleMode();
  void setCreateMode();
  void setEditMode();

  short last_x, last_y;
  short curr_x, curr_y;
  MEvent::MouseButtonType mouseButton;
  MEvent::ModifierType modifier;

  std::shared_ptr<apps::interaction::InteractionMode> interactionMode;
  std::shared_ptr<canvas::SquidgetToolManager> squidgetToolManager;
  std::shared_ptr<apps::CurveHUD> curveHUD;

  std::shared_ptr<common::logging::Logger> logger;
};
} // namespace apps
} // namespace bookmark_squidget

#endif