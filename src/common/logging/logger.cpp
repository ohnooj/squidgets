#include "common/logging/logger.hpp"

namespace common {
namespace logging {
std::string timestep() {
  auto currentTime =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

  // Convert the current time to a tm structure for formatting
  std::tm *timeInfo = std::localtime(&currentTime);

  // Format the time and store it in a string
  std::stringstream formattedTime;
  formattedTime << std::setw(4) << std::setfill('0')
                << (timeInfo->tm_year + 1900) << ":" // Year
                << std::setw(2) << std::setfill('0') << (timeInfo->tm_mon + 1)
                << ":" // Month
                << std::setw(2) << std::setfill('0') << timeInfo->tm_mday
                << ":" // Day
                << std::setw(2) << std::setfill('0') << timeInfo->tm_hour
                << ":" // Hour
                << std::setw(2) << std::setfill('0') << timeInfo->tm_min
                << ":" // Minute
                << std::setw(2) << std::setfill('0')
                << timeInfo->tm_sec // Second
                << ":"              // Minute
                << std::setw(3) << std::setfill('0')
                << (std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::system_clock::now().time_since_epoch()) %
                    1000)
                       .count(); // Milliseconds

  // Retrieve the formatted time as a string
  std::string formattedTimeString = formattedTime.str();
  return formattedTimeString;
}

// Setup Logs ==================================================================
void Logger::logSystem(const std::string &message) { log("SYSTEM", message); }

void Logger::logPluginLoaded() { logSystem("plugin-loaded," + loggerName); }

void Logger::logPluginActivated() { logSystem("plugin-activated," + loggerName); }

void Logger::logPluginDeactivated() { logSystem("plugin-deactivated," + loggerName); }

void Logger::logPluginUnloaded() { logSystem("plugin unloaded," + loggerName); }

// Tool ========================================================================
void Logger::logTool(const std::string &message) { log("TOOL", message); }

void Logger::logInputMode(const int inputMode) {
  // InputModes:
  // 0 - free
  // 1 - 1 Stroke
  // 2 - 1 Stroke + live
  // 3 - 2 Stroke
  std::string im = "";
  if (inputMode == 0) {
    im = "freemode";
  } else if (inputMode == 1) {
    im = "1 Stroke";
  } else if (inputMode == 2) {
    im = "1 Stroke + live";
  } else if (inputMode == 3) {
    im = "2 Stroke";
  }

  logTool("INPUT MODE," + im);
}

void Logger::logTaskMode(const int taskMode) {
  // TaskModes:
  // 0 - free
  // 1 - translate
  // 2 - T + R
  // 3 - Deform
  std::string tm = "";
  if (taskMode == 0) {
    tm = "no task";
  } else if (taskMode == 1) {
    tm = "translate";
  } else if (taskMode == 2) {
    tm = "T + R";
  } else if (taskMode == 3) {
    tm = "Deform";
  }

  logTool("TASK MODE," + tm);
}

void Logger::logSketchStroke(const std::string &stroke) {
  logTool("SKETCH," + stroke);
}

void Logger::logProjectedCurve(const std::string &point) {
  logTool("PROJECTED," + point);
}

void Logger::logSelectedObjectName(const std::string &objectName) {
  logTool("SELECTED," + objectName);
}

void Logger::logSquidgetType(const std::string &squidgetType) {
  logTool("SQUIDGET-TYPE," + squidgetType);
}

// Canvas ======================================================================
void Logger::logCanvas(const std::string &message) { log("CANVAS", message); }

// Bookmark ====================================================================
void Logger::logBookmark(const std::string &message) { log("BOOKMARK", message); }

void Logger::logCreateMode() { logBookmark("CREATE-MODE"); }

void Logger::logEditMode() { logBookmark("EDIT-MODE"); }


// Geometry ====================================================================
void Logger::logGeometry(const std::string &message) {
  log("GEOMETRY", message);
}

void Logger::logTransformMatrix(const std::string obj,
                                const MMatrix &old_matrix,
                                const MMatrix &new_matrix) {
  std::string old_matrix_str = "";
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      old_matrix_str += std::to_string(old_matrix(i, j)) + " ";
    }
  }
  old_matrix_str.pop_back();

  std::string new_matrix_str = "";
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      new_matrix_str += std::to_string(new_matrix(i, j)) + " ";
    }
  }
  new_matrix_str.pop_back();

  std::string together =
      "MATRIX," + obj + "," + old_matrix_str + "," + new_matrix_str;
  logGeometry(together);
}

void Logger::logVertexPositions(const std::string obj,
                                const MPointArray &oldVertices,
                                const MPointArray &newVertices) {
  std::string oldVertices_str = "";
  for (unsigned int i = 0; i < oldVertices.length(); i++) {
    oldVertices_str += std::to_string(oldVertices[i].x) + " " +
                       std::to_string(oldVertices[i].y) + " " +
                       std::to_string(oldVertices[i].z) + " ";
  }
  oldVertices_str.pop_back();

  std::string newVertices_str = "";
  for (unsigned int i = 0; i < newVertices.length(); i++) {
    newVertices_str += std::to_string(newVertices[i].x) + " " +
                       std::to_string(newVertices[i].y) + " " +
                       std::to_string(newVertices[i].z) + " ";
  }
  newVertices_str.pop_back();

  std::string together =
      "VERTEX," + obj + "," + oldVertices_str + "," + newVertices_str;
  logGeometry(together);
}

void Logger::logCandidateOldShadow(const std::string obj,
                                   const MPointArray &shadowWorldCVs) {
  std::string shadowWorldCVs_str = "";
  for (unsigned int i = 0; i < shadowWorldCVs.length(); i++) {
    shadowWorldCVs_str += std::to_string(shadowWorldCVs[i].x) + " " +
                          std::to_string(shadowWorldCVs[i].y) + " " +
                          std::to_string(shadowWorldCVs[i].z) + " ";
  }
  shadowWorldCVs_str.pop_back(); // remove last space

  std::string together = "CANDIDATE," + obj + "," + shadowWorldCVs_str;
  logGeometry(together);
}

void Logger::logCurveDistScore(const std::string obj, const double score,
                               const double squared_norm,
                               const double matrix_distance) {
  std::string together = "SCORE," + obj + "," + std::to_string(score) + "," +
                         std::to_string(squared_norm) + "," +
                         std::to_string(matrix_distance);
  logGeometry(together);
}

// Command
// =====================================================================

void Logger::logCommand(const std::string &command) { log("COMMAND", command); }

// Mouse IO
// ====================================================================
void Logger::logMouse(const std::string &message) { log("MOUSE", message); }

void Logger::logMousePress(int mouseButton, int x, int y) {
  std::string msg = "";
  if (mouseButton == 0) {
    msg = "left";
  } else if (mouseButton == 1) {
    msg = "middle";
  }
  msg += "," + std::to_string(x) + " " + std::to_string(y);

  logMouse("press," + msg);
}

void Logger::logMouseHold() { logMouse("hold"); }

void Logger::logMouseTap() { logMouse("tap"); }

void Logger::logMouseRelease(int mouseButton, int x, int y) {
  std::string msg = "";
  if (mouseButton == 0) {
    msg = "left";
  } else if (mouseButton == 1) {
    msg = "middle";
  }
  msg += "," + std::to_string(x) + " " + std::to_string(y);

  logMouse("release");
}

// Keyboard IO
// =================================================================
void Logger::logKeyboard(const std::string &message) { log("KEY", message); }

void Logger::logAbortAction() { logKeyboard("abort"); }

void Logger::logCompleteAction() { logKeyboard("complete"); }

void Logger::logDeleteAction() { logKeyboard("delete"); }

void Logger::logUndo() { logKeyboard("undo"); }

void Logger::logShift() { logKeyboard("shift"); }

void Logger::logCtrl() { logKeyboard("ctrl"); }

// Logging
// =====================================================================
void Logger::log(const std::string &type, const std::string &message) {
  // cout << timestep() << "," << type << "," << message << endl;

  // Calling the abstractSquidgetToolCmd raw breaks the logger for some reason.
  // so just recompile the tool and return.  I don't want to deal with this now.
  // return

  if (logFilepath.empty()) {
    // skip logging
    cout << "No log file path set" << endl;
    return;
  }

  std::ofstream outfile(logFilepath, std::ios::app);

  if (!outfile.is_open()) {
    cout << "THROWING ERROR" << endl;
    std::cerr << "Failed to open log file: " << logFilepath << std::endl;
  }
  outfile << timestep() << "," << type << "," << message << std::endl;
  outfile.close();
}

} // namespace logging
} // namespace common