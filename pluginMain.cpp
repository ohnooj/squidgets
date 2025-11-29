#include <fstream>
#include <iostream>

#include <maya/MFnPlugin.h>
#include <maya/MStreamUtils.h>

#include "abstract_squidget/apps/abstractContextCmd.hpp"
#include "abstract_squidget/apps/abstractSquidgetTool.hpp"
#include "bookmark_squidget/apps/bookmarkContextCmd.hpp"
#include "bookmark_squidget/apps/bookmarkSquidgetTool.hpp"

#define ABSTRACT_CONTEXT_NAME "ac316"
#define ABSTRACT_TOOL_NAME "at316"
#define ABSTRACT_CONTEXT "abstractContext"
#define ABSTRACT_TOOL_CMD "abstractSquidgetToolCmd"

#define BOOKMARK_CONTEXT_NAME "bc316"
#define BOOKMARK_TOOL_NAME "bt316"
#define BOOKMARK_CONTEXT "bookmarkContext"
#define BOOKMARK_TOOL_CMD "bookmarkSquidgetToolCmd"

void createButtons();
void deleteButtons();

using namespace abstract_squidget::apps;
using namespace bookmark_squidget::apps;

MStatus initializePlugin(MObject obj) {
  // Cout to standard out
  cout.rdbuf(MStreamUtils::stdOutStream().rdbuf());
  std::cout.rdbuf(MStreamUtils::stdOutStream().rdbuf());
  MStatus status;

  // Squidget Plugin
  MFnPlugin aplugin(obj, "Abstract Squidget", "0.1", "Any");
  status = aplugin.registerContextCommand(
      ABSTRACT_CONTEXT, AbstractContextCmd::creator, ABSTRACT_TOOL_CMD,
      AbstractSquidgetTool::creator, AbstractSquidgetTool::newSyntax);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  // Bookmark Creation Plugin
  MFnPlugin bplugin(obj, "Bookmark Squidget", "0.1", "Any");
  status = bplugin.registerContextCommand(
      BOOKMARK_CONTEXT, BookmarkContextCmd::creator, BOOKMARK_TOOL_CMD,
      BookmarkSquidgetTool::creator, BookmarkSquidgetTool::newSyntax);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  createButtons();

  return status;
}

MStatus uninitializePlugin(MObject obj) {
  deleteButtons();
  MStatus status;
  MFnPlugin plugin(obj);

  status = plugin.deregisterContextCommand(ABSTRACT_CONTEXT, ABSTRACT_TOOL_CMD);
  status = plugin.deregisterContextCommand(BOOKMARK_CONTEXT, BOOKMARK_TOOL_CMD);

  CHECK_MSTATUS_AND_RETURN_IT(status);
  return status;
}

void createButtons() {
  std::ostringstream ss;
  ss << ABSTRACT_CONTEXT " " ABSTRACT_CONTEXT_NAME << ";" << endl;
  ss << "setParent Custom;" << endl;
  ss << "toolButton" << endl;
  ss << "    -cl toolCluster" << endl;
  ss << "    -enableBackground true" << endl;
  ss << "    -bgc 0.1039 0.2353 0.3843" << endl;
  ss << "    -t " ABSTRACT_CONTEXT_NAME << endl;
  ss << "    -i1 \"commandButton.png\" " << endl;
  ss << "    " ABSTRACT_TOOL_NAME << ";" << endl;
  MGlobal::executeCommand(ss.str().c_str());

  ss.str("");
  ss << BOOKMARK_CONTEXT " " BOOKMARK_CONTEXT_NAME << ";" << endl;
  ss << "setParent Custom;" << endl;
  ss << "toolButton" << endl;
  ss << "    -cl toolCluster" << endl;
  ss << "    -enableBackground true" << endl;
  ss << "    -bgc 0.6706 0.3451 0.7765" << endl;
  ss << "    -t " BOOKMARK_CONTEXT_NAME << endl;
  ss << "    -i1 \"commandButton.png\" " << endl;
  ss << "    " BOOKMARK_TOOL_NAME << ";" << endl;
  MGlobal::executeCommand(ss.str().c_str());
}

void deleteButtons() {
  MGlobal::executeCommand("flushUndo; deleteUI " ABSTRACT_CONTEXT_NAME "; "
                          "deleteUI " ABSTRACT_TOOL_NAME "; flushUndo;");
  MGlobal::executeCommand("flushUndo; deleteUI " BOOKMARK_CONTEXT_NAME "; "
                          "deleteUI " BOOKMARK_TOOL_NAME "; flushUndo;");
}
