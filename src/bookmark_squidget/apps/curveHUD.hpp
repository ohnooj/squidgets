#ifndef CURVEHUD_HPP
#define CURVEHUD_HPP

#include <iostream>
#include <sstream>

#include <maya/MGlobal.h>
#include <maya/MString.h>

namespace bookmark_squidget {
namespace apps {
class CurveHUD {
public:
  CurveHUD();
  ~CurveHUD() {}

  void drawHUD();
  void removeHUD();
  void updateHUD();

  void setMode(bool mode) { isCreateMode = mode; }
  void setCanvasVisible(bool visible) { isEditCanvasVisible = visible; }

  void drawCreateDrewOnCanvas(bool drawn);

private:
  bool isCreateMode;
  bool drewOnCanvas;
  bool isEditCanvasVisible;

  MString modeHUD;
  MString canvasHUD;
  // MString listHUD;

  void setupHUD();
  void drawCurveMode();
  void drawCanvasVisible();
  // void drawSelectionList();
  void showGrid(bool show);
  void showControls();

  void updateModeHUD();
  void updateCanvasVisibleHUD();
  // void updateListHUD();

  // void removeToolListHUD();
  void removeEditCanvasVisibilityHUD();
};
} // namespace apps
} // namespace bookmark_squidget

#endif // CURVEHUD_HPP