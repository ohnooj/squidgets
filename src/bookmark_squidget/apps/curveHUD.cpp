#include "bookmark_squidget/apps/curveHUD.hpp"

namespace bookmark_squidget {
namespace apps {

CurveHUD::CurveHUD() {
  // cout << "CurveHUD::CurveHUD" << endl;
  setupHUD();
}

// Setup HUD ================================================================
void CurveHUD::setupHUD() {
  // cout << "CurveHUD::setupHUD" << endl;
  MGlobal::executePythonCommand("import maya.cmds as cmds");
  isCreateMode = true;
  drewOnCanvas = false;
  isEditCanvasVisible = true;

  modeHUD = MString("HUDBookmarkSquidgetToolMode");
  canvasHUD = MString("HUDBookmarkSquidgetToolEditCanvasVisibility");
  // listHUD = MString("HUDBookmarkSquidgetToolListLabel");

  showControls();
  showGrid(true);
}

void CurveHUD::showControls() {
  std::ostringstream ss;
  MString modeLabel = (isCreateMode) ? "'Toggle Mode <ESC>: CREATE'"
                                     : "'Toggle Mode <ESC>: CONTROL'";
  ss << "cmds.headsUpDisplay( '" << modeHUD << "', section=1, block=0," << endl;
  ss << " blockSize='large', blockAlignment='left', " << endl;
  ss << " label= " << modeLabel << ", " << endl;
  ss << " labelWidth=300, " << endl;
  ss << " labelFontSize='large')" << endl;

  std::string cmd = ss.str();
  MGlobal::executePythonCommand(cmd.c_str());
}

void CurveHUD::showGrid(bool show) {
  // cout << "CurveHUD::showGrid" << endl;
  std::string showStr = (show) ? "True" : "False";

  std::ostringstream ss;
  ss << "cmds.headsUpDisplay( '" << modeHUD << "', edit=True, " << endl;
  ss << " showGrid=" << showStr << ")" << endl;

  std::string cmd = ss.str();
  MGlobal::executePythonCommand(cmd.c_str());
}

// Update HUD ================================================================
void CurveHUD::drawHUD() {
  // cout << "Create Mode HUD: " << isCreateMode << endl;

  drawCurveMode();
  drawCanvasVisible();
  // if (isCreateMode) {
  //   drawSelectionList();
  // }
}

void CurveHUD::drawCurveMode() {
  std::ostringstream ss;
  MString modeLabel = (isCreateMode) ? "'Toggle Mode <ESC>: CREATE'"
                                     : "'Toggle Mode <ESC>: CONTROL'";
  ss << "cmds.headsUpDisplay( '" << modeHUD << "', section=1, block=0," << endl;
  ss << " blockSize='large', blockAlignment='left', " << endl;
  ss << " label= " << modeLabel << ", " << endl;
  ss << " labelWidth=300, " << endl;
  ss << " labelFontSize='large')" << endl;

  std::string cmd = ss.str();
  // cout << "Command: " << cmd << endl;
  MGlobal::executePythonCommand(cmd.c_str());
  updateModeHUD();
}

// void CurveHUD::drawSelectionList() {
//   if (!isCreateMode) return;
//   std::ostringstream ss;
//   ss << "def objectPosition(*args):" << endl;
//   // ss << "    global BookmarkSquidgetToolSelNodes" << endl;
//   ss << "    try:" << endl;
//   ss << "        i = 0" << endl;
//   ss << "        row_name = 'HUDBookmarkSquidgetToolListItem' + str(i)" <<
//   endl; ss << "        HUDList = cmds.headsUpDisplay(lh=True)" << endl; ss <<
//   "        while row_name in HUDList:" << endl; ss << "
//   cmds.headsUpDisplay(row_name, remove=True)" << endl; ss << "            i =
//   i + 1" << endl; ss << "            row_name =
//   'HUDBookmarkSquidgetToolListItem' + str(i)" << endl; ss << "
//   BookmarkSquidgetToolSelNodes = cmds.selectedNodes(dagObjects=True)" <<
//   endl; ss << "        BookmarkSquidgetToolSelNodes = [n.split('|')[-1] for n
//   in BookmarkSquidgetToolSelNodes]" << endl; ss << "        for i in
//   range(len(BookmarkSquidgetToolSelNodes)):" << endl; ss << "            node
//   = BookmarkSquidgetToolSelNodes[i]" << endl; ss << "            row_name =
//   'HUDBookmarkSquidgetToolListItem' + str(i)" << endl; ss << "
//   cmds.headsUpDisplay( row_name, section=1, block=3+i, " << endl; ss << "
//   blockSize='small', blockAlignment='left', label=node, " << endl; ss << "
//   labelFontSize='small')" << endl; ss << "        return ''" << endl; ss << "
//   except:" << endl; ss << "        return ''" << endl; ss <<
//   "cmds.headsUpDisplay( '" << listHUD << "', section=1, block=2, " << endl;
//   ss << "    blockSize='small', blockAlignment='left', label='Selected
//   Objects to Key: click (+ <shift>)', " << endl; ss << "    labelWidth=300, "
//   << endl; ss << "    labelFontSize='small', command=objectPosition,
//   event='SelectionChanged'," << endl; ss << " nodeChanges='attributeChange'
//   )" << endl; std::string cmd = ss.str();
//   // cout << "Command: " << cmd << endl;
//   MGlobal::executePythonCommand(cmd.c_str());
// }

void CurveHUD::drawCanvasVisible() {
  std::ostringstream ss;
  ss << "cmds.headsUpDisplay( '" << canvasHUD << "', section=1, block=1, "
     << endl;
  ss << "    blockSize='medium', blockAlignment='left', label='Toggle Canvas "
        "Visibility <enter>', "
     << endl;
  ss << "    labelWidth=300, " << endl;
  ss << "    labelFontSize='large' )" << endl;

  std::string cmd = ss.str();
  // cout << "Command: " << cmd << endl;
  MGlobal::executePythonCommand(cmd.c_str());

  updateCanvasVisibleHUD();
}

// Update HUD ================================================================
void CurveHUD::updateHUD() {
  updateModeHUD();
  updateCanvasVisibleHUD();
  // updateListHUD();
}

void CurveHUD::updateModeHUD() {
  std::ostringstream ss;
  MString modeLabel = (isCreateMode) ? "'Toggle Mode <ESC>: CREATE'"
                                     : "'Toggle Mode <ESC>: CONTROL'";

  ss << "cmds.headsUpDisplay( '" << modeHUD << "', edit=True, " << endl;
  ss << " label= " << modeLabel << ")" << endl;

  std::string cmd = ss.str();
  // cout << "Command: " << cmd << endl;
  MGlobal::executePythonCommand(cmd.c_str());
}

void CurveHUD::updateCanvasVisibleHUD() {
  std::ostringstream ss;
  MString canvasVisible = (isEditCanvasVisible) ? "" : "HIDDEN";
  MString canvasLabel = (isCreateMode)
                            ? "'Create Canvas for Selected Objects <enter>'"
                            : "'Canvas Visibility <enter>: " + canvasVisible + "'";

  ss << "cmds.headsUpDisplay( '" << canvasHUD << "', edit=True, " << endl;
  ss << " label= " << canvasLabel << ")" << endl;

  std::string cmd = ss.str();
  // cout << "Command: " << cmd << endl;
  MGlobal::executePythonCommand(cmd.c_str());
}

// void CurveHUD::updateListHUD() {
//   MString visible = (isCreateMode) ? "True" : "False";
//   std::ostringstream ss;
//   ss << "i = 0" << endl;
//   ss << "cmds.headsUpDisplay( '" << listHUD << "', edit=True, " << endl;
//   ss << "    visible= " << visible << ")" << endl;
//   ss << "row_name = 'HUDBookmarkSquidgetToolListItem' + str(i)" << endl;
//   ss << "HUDList = cmds.headsUpDisplay(lh=True)" << endl;
//   ss << "while row_name in HUDList:" << endl;
//   ss << "    cmds.headsUpDisplay(row_name, edit=True, " << endl;
//   ss << "        visible=" << visible <<  " )" << endl;
//   ss << "    i = i + 1" << endl;
//   ss << "    row_name = 'HUDBookmarkSquidgetToolListItem' + str(i)" << endl;

//   std::string cmd = ss.str();
//   // cout << "Command: " << cmd << endl;
//   MGlobal::executePythonCommand(cmd.c_str());
// }

// Remove HUD ================================================================
void CurveHUD::removeHUD() {
  MGlobal::executePythonCommand("cmds.headsUpDisplay('" + modeHUD +
                                "', remove=True)");
  removeEditCanvasVisibilityHUD();
  // removeToolListHUD();
}

void CurveHUD::removeEditCanvasVisibilityHUD() {
  MGlobal::executePythonCommand("cmds.headsUpDisplay('" + canvasHUD +
                                "', remove=True)");
}

// void CurveHUD::removeToolListHUD() {
//   if (!isCreateMode) return;

//   MGlobal::executePythonCommand("cmds.headsUpDisplay('" + listHUD +  "',
//   remove=True)");

//   std::ostringstream ss;
//   ss << "i = 0" << endl;
//   ss << "row_name = 'HUDBookmarkSquidgetToolListItem' + str(i)" << endl;
//   ss << "HUDList = cmds.headsUpDisplay(lh=True)" << endl;
//   ss << "while row_name in HUDList:" << endl;
//   ss << "    cmds.headsUpDisplay(row_name, remove=True)" << endl;
//   ss << "    i = i + 1" << endl;
//   ss << "    row_name = 'HUDBookmarkSquidgetToolListItem' + str(i)" << endl;

//   std::string cmd = ss.str();
//   // cout << "Command: " << cmd << endl;
//   MGlobal::executePythonCommand(cmd.c_str());
// }

void CurveHUD::drawCreateDrewOnCanvas(bool drawn) {
  drewOnCanvas = drawn;
  if (!drewOnCanvas) {
    std::ostringstream ss;
    ss << "cmds.headsUpDisplay( 'HUDBookmarkSquidgetToolCanvasDraw', "
          "section=1, block=3,"
       << endl;
    ss << " blockSize='large', blockAlignment='left', " << endl;
    ss << " label='Stroke was not drawn on a canvas', " << endl;
    ss << " labelFontSize='large')" << endl;

    std::string cmd = ss.str();
    // cout << "Command: " << cmd << endl;
    MGlobal::executePythonCommand(cmd.c_str());
  } else {
    MGlobal::executePythonCommand(
        "cmds.headsUpDisplay('HUDBookmarkSquidgetToolCanvasDraw', "
        "remove=True)");
  }
}

} // namespace apps
} // namespace bookmark_squidget