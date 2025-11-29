/**
 * @file logger.hpp
 * @brief This file defines the Logger class.
 * We use the Logger class to log user input into a file.
 * The file created is named: "squidget-yyyymmdd-hhmmss.log"
 *
 * The format of the file is as follows:
 * - [yyyy:mm:dd:hh:mm:ss]: <message>
 *
 * Where message can be any of the following:
 * - "plugin-loaded"
 * - "plugin-activated"
 * - "plugin-deactivated"
 *
 * Create mode
 * - "canv-cre: <canvas_id>"
 * - "canv-del: <canvas_id>"
 * - "squi-cre: <canvas_id>:<squidget_id>"
 * - "squi-del: <canvas_id>:<squidget_id>"
 * - "squi-con: <canvas_id>:<squidget_id>:<squidget_id>:<squidget_id>..."
 * - "squi-dis: <canvas_id>:<squidget_id>:<squidget_id>:<squidget_id>..."
 *
 * - "SQUIDGET_MODE: <mode>"
 * - "disc-set: <canvas_id>:<squidget_id>:<value>"
 * - "cont-set: <canvas_id>:<squidget_id>:<value>"
 * - "cont-hol: <canvas_id>:<squidget_id>:<value>"
 * - "cont-rel: <canvas_id>:<squidget_id>:<value>"
 *
 * - "impl-set: <name>"
 * - "impl-hol: <name>:<value>"
 * - "impl-rel: <name>:<value>"
 * - "plugin unloaded"
 *
 * @date 2023-11-18
 *
 * @copyright Copyright (c) 2023
 */

#ifndef LOGGER_HPP_INCLUDED
#define LOGGER_HPP_INCLUDED

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <maya/MArgList.h>
#include <maya/MMatrix.h>
#include <maya/MPointArray.h>

namespace common {
namespace logging {

std::string timestep();

class Logger {
public:
  Logger() {
    logFilepath = "";
    loggerName = "No-Logger-Name";
    cout << "Logger():<" << logFilepath << ", " << loggerName << ">" << endl;
  };

  Logger(std::string logFileName, std::string loggerName) {
    logFilepath = logFileName;
    loggerName = loggerName;
    cout << "Logger():<" << logFilepath << ", " << loggerName << ">" << endl;
  }

  ~Logger() {};

  // System ========================================================
  void logPluginLoaded();
  void logPluginUnloaded();
  void logPluginActivated();   // ./
  void logPluginDeactivated(); // ./

  // Tool ==========================================================
  void logInputMode(const int inputMode);                    // ./
  void logTaskMode(const int taskMode);                      // ./
  void logSketchStroke(const std::string &stroke);           // ./
  void logProjectedCurve(const std::string &point);          // ./
  void logSelectedObjectName(const std::string &objectName); // ./
  void logSquidgetType(const std::string &squidgetType);     // ./

  // Bookmark ======================================================
  void logCreateMode(); // ./
  void logEditMode(); // ./

  // Geometry ======================================================
  void logTransformMatrix(const std::string obj, const MMatrix &old_matrix,
                          const MMatrix &new_matrix); // ./
  void logVertexPositions(const std::string obj, const MPointArray &oldVertices,
                          const MPointArray &newVertices); // ./
  void logCandidateOldShadow(const std::string obj,
                             const MPointArray &shadowWorldCVs); // ./
  void logCurveDistScore(const std::string obj, const double score,
                         const double squared_norm,
                         const double matrix_distance); // ./

  // CommandTool ===================================================
  void logCommand(const std::string &command); // ./

  // Mouse IO ======================================================
  void logMousePress(int mouseButton, int x, int y);   // ./
  void logMouseHold();                                 // ./
  void logMouseTap();                                  // ./
  void logMouseRelease(int mosueButton, int x, int y); // ./

  // Keyboard IO ===================================================
  void logCompleteAction(); // ./
  void logDeleteAction(); // ./
  void logAbortAction(); // ./
  void logUndo();        // ./

  void logShift(); // ./
  void logCtrl();  // ./

  void setLogFilepath(const std::string &logFilepath) {
    cout << "Logger::setLogFilepath" << logFilepath << endl;
    this->logFilepath = logFilepath;
  };

  std::string getLogFilepath() { return logFilepath; };

  void log(std::string message) {
    cout << "Logger::log(): " << message << endl;
    logSystem(message);
  };

private:
  std::string logFilepath;
  std::string loggerName;

  void logSystem(const std::string &message);
  void logTool(const std::string &message);
  void logBookmark(const std::string &message);
  void logCanvas(const std::string &message);
  void logGeometry(const std::string &message);
  void logMouse(const std::string &button);
  void logKeyboard(const std::string &message);
  void log(const std::string &type, const std::string &message);
};
} // namespace logging
} // namespace common
#endif // LOGGER_HPP_INCLUDED
