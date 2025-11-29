#include "abstract_squidget/apps/abstractContextCmd.hpp"

#include <maya/MArgList.h>
#include <maya/MStringArray.h>

namespace abstract_squidget {
namespace apps {
AbstractContextCmd::AbstractContextCmd() {
  cout << "AbstractContextCmd()" << endl;
}

AbstractContextCmd::~AbstractContextCmd() {
  cout << "~AbstractContextCmd()" << endl;
  delete fAbstractContext;
}

MStatus AbstractContextCmd::doEditFlags() {
  MStatus status = MS::kSuccess;
  MArgParser argData = parser();

  if (argData.isFlagSet(kStrokePointDistFlag)) {
    unsigned strokePointDist;
    status = argData.getFlagArgument(kStrokePointDistFlag, 0, strokePointDist);
    if (!status) {
      status.perror("minstrokePointDist flag parsing failed.");
      return status;
    }
    fAbstractContext->setStrokePointDist(strokePointDist);
  }
  if (argData.isFlagSet(kStrokeWidthFlag)) {
    unsigned strokeWidth;
    status = argData.getFlagArgument(kStrokeWidthFlag, 0, strokeWidth);
    if (!status) {
      status.perror("strokeWidth flag parsing failed.");
      return status;
    }
    fAbstractContext->setStrokeWidth(strokeWidth);
  }
  if (argData.isFlagSet(kInputModeFlag)) {
    unsigned mode;
    status = argData.getFlagArgument(kInputModeFlag, 0, mode);
    if (!status) {
      status.perror("mode flag parsing failed.");
      return status;
    }
    fAbstractContext->setInputMode(mode);
  }
  if (argData.isFlagSet(kTaskModeFlag)) {
    unsigned mode;
    status = argData.getFlagArgument(kTaskModeFlag, 0, mode);
    if (!status) {
      status.perror("mode flag parsing failed.");
      return status;
    }
    fAbstractContext->setTaskMode(mode);
  }
  if (argData.isFlagSet(kCandidateObjectsFlag)) {
    MString objCandidates;
    status = argData.getFlagArgument(kCandidateObjectsFlag, 0, objCandidates);
    if (!status) {
      status.perror("objCandidates flag parsing failed.");
      return status;
    }
    fAbstractContext->setObjCandidatesStr(objCandidates);
  }
  if (argData.isFlagSet(kLoggerPathFlag)) {
    MString loggerPath;
    status = argData.getFlagArgument(kLoggerPathFlag, 0, loggerPath);
    if (!status) {
      status.perror("loggerPath flag parsing failed.");
      return status;
    }
    fAbstractContext->setLoggerPath(loggerPath);
  }
  if (argData.isFlagSet(kSquidgetTypeFlag)) {
    MString squidgetType;
    status = argData.getFlagArgument(kSquidgetTypeFlag, 0, squidgetType);
    if (!status) {
      status.perror("squidgetType flag parsing failed.");
      return status;
    }
    fAbstractContext->setSquidgetType(squidgetType);
  }
  return status;
}

MStatus AbstractContextCmd::doQueryFlags() {
  MArgParser argData = parser();
  if (argData.isFlagSet(kStrokePointDistFlag)) {
    setResult((int)fAbstractContext->getStrokePointDist());
  }
  if (argData.isFlagSet(kStrokeWidthFlag)) {
    setResult((int)fAbstractContext->getStrokeWidth());
  }
  if (argData.isFlagSet(kInputModeFlag)) {
    setResult((int)fAbstractContext->getInputMode());
  }
  if (argData.isFlagSet(kTaskModeFlag)) {
    setResult((int)fAbstractContext->getTaskMode());
  }
  if (argData.isFlagSet(kCandidateObjectsFlag)) {
    setResult((MString)fAbstractContext->getObjCandidatesStr());
  }
  if (argData.isFlagSet(kLoggerPathFlag)) {
    setResult((MString)fAbstractContext->getLoggerPath());
  }
  if (argData.isFlagSet(kLoggerPathFlag)) {
    setResult((MString)fAbstractContext->getLoggerPath());
  }
  if (argData.isFlagSet(kSquidgetTypeFlag)) {
    setResult((MString)fAbstractContext->getSquidgetType());
  }
  return MS::kSuccess;
}

MPxContext *AbstractContextCmd::makeObj() {
  // cout << "AbstractContextCmd::makeObj()" << endl;
  fAbstractContext = new AbstractContext();
  return fAbstractContext;
}

MStatus AbstractContextCmd::appendSyntax() {
  MSyntax mySyntax = syntax();
  if (MS::kSuccess != mySyntax.addFlag(kStrokePointDistFlag,
                                       kStrokePointDistFlagLong,
                                       MSyntax::kUnsigned)) {
    return MS::kFailure;
  }
  if (MS::kSuccess != mySyntax.addFlag(kStrokeWidthFlag, kStrokeWidthFlagLong,
                                       MSyntax::kUnsigned)) {
    return MS::kFailure;
  }
  if (MS::kSuccess != mySyntax.addFlag(kInputModeFlag, kInputModeFlagLong,
                                       MSyntax::kUnsigned)) {
    return MS::kFailure;
  }
  if (MS::kSuccess !=
      mySyntax.addFlag(kTaskModeFlag, kTaskModeFlagLong, MSyntax::kUnsigned)) {
    return MS::kFailure;
  }
  if (MS::kSuccess != mySyntax.addFlag(kCandidateObjectsFlag,
                                       kCandidateObjectsFlagLong,
                                       MSyntax::kString)) {
    return MS::kFailure;
  }
  if (MS::kSuccess != mySyntax.addFlag(kLoggerPathFlag, kLoggerPathFlagLong,
                                       MSyntax::kString)) {
    return MS::kFailure;
  }
  if (MS::kSuccess != mySyntax.addFlag(kSquidgetTypeFlag, kSquidgetTypeFlagLong,
                                       MSyntax::kString)) {
    return MS::kFailure;
  }

  return MS::kSuccess;
}

void *AbstractContextCmd::creator() { return new AbstractContextCmd; }

} // namespace apps
} // namespace abstract_squidget