#include "abstract_squidget/apps/abstractSquidgetTool.hpp"

#include <queue>
#include <random>

#include <maya/MArgDatabase.h>
#include <maya/MArgList.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MSyntax.h>

#include "abstract_squidget/deformer/abstractDeformer.hpp"
#include "abstract_squidget/deformer/optimization.hpp"
#include "common/maya/toonRenderer.hpp"
#include "common/util/math.hpp"
#include "common/util/maya.hpp"

namespace abstract_squidget {
namespace apps {
AbstractSquidgetTool::AbstractSquidgetTool() {
  setCommandString("abstractSquidgetTool");
  setHistoryOn(true);

  setPreSelect(false);
  setDeform(false);
}

AbstractSquidgetTool::~AbstractSquidgetTool() {}

// =============================================================================
//                      Tool Command Typical Methods
// =============================================================================
MSyntax AbstractSquidgetTool::newSyntax() {
  // cout << "AbstractSquidgetTool::newSyntax()" << endl;
  MSyntax syntax;
  syntax.addFlag(kCurveFlag, kCurveFlagLong, MSyntax::kString);
  syntax.addFlag(kCandidateObjectsFlag, kCandidateObjectsFlagLong,
                 MSyntax::kString);
  syntax.addFlag(kPreSelectFlag, kPreSelectFlagLong, MSyntax::kBoolean);
  syntax.addFlag(kDeformFlag, kDeformFlagLong, MSyntax::kBoolean);
  syntax.addFlag(kTranslateOnlyFlag, kTranslateOnlyFlagLong, MSyntax::kBoolean);
  syntax.addFlag(kShadowCurveFlag, kShadowCurveFlagLong, MSyntax::kString);
  return syntax;
}

void *AbstractSquidgetTool::creator() {
  // cout << "AbstractSquidgetTool::creator()" << endl;
  return new AbstractSquidgetTool;
}

MStatus AbstractSquidgetTool::doIt(const MArgList &args) {
  // cout << "AbstractSquidgetTool::doIt()" << endl;
  MStatus status = parseArgs(args);

  if (MS::kSuccess != status) {
    return status;
  }
  return redoIt();
}

MStatus AbstractSquidgetTool::undoIt() {
  MStatus status = MS::kSuccess;

  MDagPath geomPath = common::util::getDagPathFromName(editedObject);
  if (!preSelect) { // Preselct has no manipulation
    if (deform) {
      MFnMesh fnMesh(geomPath);
      fnMesh.setPoints(oldVertices, MSpace::kObject);
    } else {
      MFnTransform fnTransform(geomPath, &status);
      fnTransform.set(oldMatrix);
    }
  }
  // Set old editedObjStr and shadowCurvesStr
  *editedObj = oldEditedObj;
  *shadowCurveLocalStr = oldShadowCurveLocalStr;
  *shadowCurveWorld = oldShadowCurveWorld;

  logger->logUndo();

  return status;
}

MStatus AbstractSquidgetTool::parseArgs(const MArgList &args) {
  // cout << "AbstractSquidgetTool::parseArgs()" << endl;
  MStatus status;
  MArgDatabase argData(syntax(), args);

  if (argData.isFlagSet(kCurveFlag)) {
    MString curvePoints2D;
    status = argData.getFlagArgument(kCurveFlag, 0, curvePoints2D);
    if (MS::kSuccess != status) {
      return status;
    }
    setCurvePoints2D(curvePoints2D);
  }
  if (argData.isFlagSet(kCandidateObjectsFlag)) {
    MString candidateObjs;
    status = argData.getFlagArgument(kCandidateObjectsFlag, 0, candidateObjs);
    if (MS::kSuccess != status) {
      return status;
    }
    setCandidateObjs(candidateObjs);
  }
  if (argData.isFlagSet(kPreSelectFlag)) {
    bool preSelect;
    status = argData.getFlagArgument(kPreSelectFlag, 0, preSelect);
    if (MS::kSuccess != status) {
      return status;
    }
    setPreSelect(preSelect);
  }
  if (argData.isFlagSet(kDeformFlag)) {
    bool deform;
    status = argData.getFlagArgument(kDeformFlag, 0, deform);
    if (MS::kSuccess != status) {
      return status;
    }
    setDeform(deform);
  }
  if (argData.isFlagSet(kTranslateOnlyFlag)) {
    bool translateOnly;
    status = argData.getFlagArgument(kTranslateOnlyFlag, 0, translateOnly);
    if (MS::kSuccess != status) {
      return status;
    }
    setTranslateOnly(translateOnly);
  }
  if (argData.isFlagSet(kShadowCurveFlag)) {
    MString shadowCurve;
    status = argData.getFlagArgument(kShadowCurveFlag, 0, shadowCurve);
    if (MS::kSuccess != status) {
      return status;
    }

    setCurveToMatchStr(shadowCurve);
  }

  return status;
}

MStatus AbstractSquidgetTool::cancel() {
  // cout << "AbstractSquidgetTool::cancel()" << endl;
  return MStatus::kSuccess;
}

bool AbstractSquidgetTool::isUndoable() const {
  return true;
}

bool AbstractSquidgetTool::hasSyntax() const {
  // cout << "AbstractSquidgetTool::hasSyntax()" << endl;
  return true;
}

MStatus AbstractSquidgetTool::finalize() {
  // cout << "AbstractSquidgetTool::finalize()" << endl;
  MArgList command = getFinalCommand();

  printCommand(command);
  return MPxToolCommand::doFinalize(command);
}

MArgList AbstractSquidgetTool::getFinalCommand() {
  MArgList command;
  command.addArg(commandString());
  command.addArg(MString(kCurveFlag));
  command.addArg(curvePoints2DStr);
  command.addArg(MString(kCandidateObjectsFlag));
  command.addArg(candidateObjsStr);
  command.addArg(MString(kPreSelectFlag));
  command.addArg(preSelect);
  command.addArg(MString(kDeformFlag));
  command.addArg(deform);
  command.addArg(MString(kTranslateOnlyFlag));
  command.addArg(translateOnly);
  command.addArg(MString(kShadowCurveFlag));
  command.addArg(inCurveToMatchStr);
  return command;
}

MString AbstractSquidgetTool::getCallableCommandStr() {
  std::ostringstream ss;
  ss << commandString().asChar() << " ";
  ss << kCurveFlag << " \"" << curvePoints2DStr.asChar() << "\" ";
  ss << kCandidateObjectsFlag << " \"" << candidateObjsStr.asChar() << "\" ";
  ss << kPreSelectFlag << " " << preSelect << " ";
  ss << kDeformFlag << " " << deform << " ";
  ss << kTranslateOnlyFlag << " " << translateOnly << " ";
  ss << kShadowCurveFlag << " \"" << inCurveToMatchStr.asChar() << "\"";

  return ss.str().c_str();
}

EditFlags AbstractSquidgetTool::getEditFlags() {
  EditFlags flags;
  flags.deform = deform;
  flags.preSelect = preSelect;
  flags.translateOnly = translateOnly;
  flags.isRigCurves = isRigCurves;
  return flags;
}

// =============================================================================
//                      Tool Command: REDOIT Main
// =============================================================================
MStatus AbstractSquidgetTool::redoIt() {
  cout << "AbstractSquidgetTool::redoIt()" << endl;
  MStatus status = MS::kSuccess;
  EditFlags editFlags = getEditFlags();
  RenderMatrices renderMats;
  RenderMatrices initRenderMats;
  QueryClosestObjectResult result;

  // We have preselected the object
  if (inCurveToMatch.length() > 0) {
    editFlags.preSelect = true;
    result = deformer::queryPreselectObject(candidateObjsStr, curvePoints2D,
                                            inCurveToMatch, editFlags);

  } else {
    result = deformer::queryClosestObject(candidateObjs, curvePoints2D,
                                          editFlags, logger);

    if (result.object == "") {
      cout << "AbstractSquidgetTool::redoIt() No object found" << endl;
      return MS::kFailure;
    }
  }

  // Push changes back to context
  RenderMatrices forWorld =
      (inCurveToMatch.length() > 0) ? result.renderMats : result.initRenderMats;
  editedObject = result.object;
  exportEditedObjectToContext(editedObject);
  exportShadowCurveLocalToContext(result.preShadowCurve);
  exportShadowCurveWorldToContext(result.preShadowCurve, forWorld);
  renderMats = result.renderMats;
  initRenderMats = result.initRenderMats;

  if (preSelect) {
    MPxToolCommand::setResult(editedObject);
    return status;
  }

  // Have object and renderMats to perform deformation
  if (deform) {
    DeformObjectResult res = deformer::deformObject(editedObject, curvePoints2D,
                                                    result.preShadowCurve,
                                                    initRenderMats, editFlags);
    oldVertices = res.oldVertices;

    exportShadowCurveLocalToContext(res.postShadowCurve);
    exportShadowCurveWorldToContext(res.postShadowCurve, initRenderMats);

    // Log the vertex positions
    MPointArray oV = res.oldVertices;
    MDagPath geomPath = common::util::getDagPathFromName(editedObject);
    MPointArray nV;
    MFnMesh fnMesh(geomPath);
    fnMesh.getPoints(nV, MSpace::kObject);

    logger->logVertexPositions(editedObject.asChar(), oV, nV);
  } else {
    DeformObjectResult res = deformer::transformObject(
        editedObject, curvePoints2D, renderMats, editFlags);

    exportShadowCurveLocalToContext(result.preShadowCurve);
    exportShadowCurveWorldToContext(result.preShadowCurve, forWorld);

    oldMatrix = res.oldTransform;

    // Log the transformation matrix
    MMatrix oM = oldMatrix.asMatrix();
    MDagPath geomPath = common::util::getDagPathFromName(editedObject);
    MFnTransform fnTransform(geomPath, &status);
    MMatrix nM = fnTransform.transformation().asMatrix();

    logger->logTransformMatrix(editedObject.asChar(), oM, nM);
  }

  MPxToolCommand::setResult(editedObject);
  return status;
}

// =============================================================================
//                      Setters for Tool Command Call
// =============================================================================
void AbstractSquidgetTool::setCurvePoints2D(MString points) {
  // Parse MString of points into MPointArray delimited by spaces
  // "0 0 0 1 1 1 2 2 2" -> MPointArray(0, 0, 0, 1, 1, 1, 2, 2, 2)
  curvePoints2DStr = points;

  MStringArray pointStrings;
  points.split(' ', pointStrings);

  // convert MStringArray to MPointArray
  curvePoints2D.clear();
  for (unsigned int i = 0; i < pointStrings.length(); i += 3) {
    MPoint p(pointStrings[i].asDouble(), pointStrings[i + 1].asDouble(),
             pointStrings[i + 2].asDouble());
    curvePoints2D.append(p);
  }
}

void AbstractSquidgetTool::setCandidateObjs(MString arr) {
  candidateObjsStr = arr;
  candidateObjs.clear();
  arr.split(' ', candidateObjs);
}

void AbstractSquidgetTool::setPreSelect(bool preSelect) {
  // cout << "AbstractSquidgetTool::setPreSelect()" << endl;
  this->preSelect = preSelect;
}

void AbstractSquidgetTool::setDeform(bool deform) {
  // cout << "AbstractSquidgetTool::setDeform()" << endl;
  this->deform = deform;
}

void AbstractSquidgetTool::setTranslateOnly(bool translateOnly) {
  // cout << "AbstractSquidgetTool::setTranslateOnly()" << endl;
  this->translateOnly = translateOnly;
}

void AbstractSquidgetTool::setCurveToMatchStr(MString shadowCurve) {
  this->inCurveToMatchStr = shadowCurve;

  MStringArray pointStrings;
  shadowCurve.split(' ', pointStrings);

  inCurveToMatch.clear();
  for (unsigned int i = 0; i < pointStrings.length(); i += 3) {
    MPoint p(pointStrings[i].asDouble(), pointStrings[i + 1].asDouble(),
             pointStrings[i + 2].asDouble());
    this->inCurveToMatch.append(p);
  }
}

void AbstractSquidgetTool::setIsRigCurves(bool rigCurves) {
  this->isRigCurves = rigCurves;
}

// =============================================================================
//                           Returns to Context
// =============================================================================
void AbstractSquidgetTool::fromContextShadowCurveWorld(
    std::shared_ptr<MPointArray> s) {
  shadowCurveWorld = s;
  oldShadowCurveWorld = *s;
}

void AbstractSquidgetTool::fromContextEditedObject(std::shared_ptr<MString> e) {
  editedObj = e;
  oldEditedObj = *e;
}

void AbstractSquidgetTool::fromContextShadowCurveLocalStr(
    std::shared_ptr<MString> s) {
  shadowCurveLocalStr = s;
  oldShadowCurveLocalStr = *s;
}

MPointArray
AbstractSquidgetTool::buildOutShadowCurveWorld(MPointArray outShadowCurve,
                                               RenderMatrices renderMats) {
  MPointArray world;
  for (MPoint p : outShadowCurve) {
    MPoint worldp = renderMats.worldMatrix * renderMats.localMatrix * p;
    world.append(worldp);
  }
  return world;
}

MString AbstractSquidgetTool::buildOutShadowCurveStr(MPointArray osc) {
  MPointArray outShadowCurve = osc;
  MString ret;
  for (unsigned int i = 0; i < outShadowCurve.length(); i++) {
    ret += outShadowCurve[i].x;
    ret += " ";
    ret += outShadowCurve[i].y;
    ret += " ";
    ret += outShadowCurve[i].z;
    ret += " ";
  }
  return ret;
}

void AbstractSquidgetTool::exportEditedObjectToContext(MString obj) {
  // cout << "AbstractSquidgetTool::getEditedObject()" << endl;
  editedObj->clear();
  *editedObj = obj;
}

void AbstractSquidgetTool::exportShadowCurveLocalToContext(MPointArray arr) {
  shadowCurveLocalStr->clear();
  *shadowCurveLocalStr = buildOutShadowCurveStr(arr);
}

void AbstractSquidgetTool::exportShadowCurveWorldToContext(MPointArray arr,
                                                           RenderMatrices rms) {
  shadowCurveWorld->clear();
  *shadowCurveWorld = buildOutShadowCurveWorld(arr, rms);
}

// =============================================================================
//                                Other
// =============================================================================
void AbstractSquidgetTool::printCommand(MArgList command) {
  MString ret;
  for (unsigned int i = 0; i < command.length(); i++) {
    ret += command.asString(i);
    ret += (i % 2 == 1) ? " " : "\n";
  }
  cout << ret << endl;
}

void AbstractSquidgetTool::saveForUndo() {
  // cout << "AbstractSquidgetTool::saveForUndo()" << endl;
  if (editedObject == "") {
    // cout << "editedObject NOT SET" << endl;
    return;
  }
  // Assumes editedObject is set
  MStatus status;
  MDagPath geomPath = common::util::getDagPathFromName(editedObject);
  MFnTransform fnTransform(geomPath, &status);
  fnTransform.set(oldMatrix);
}

void AbstractSquidgetTool::setLogger(
    std::shared_ptr<common::logging::Logger> logger) {
  // cout << "AbstractSquidgetTool::setLogger()" << endl;
  this->logger = logger;
}

} // namespace apps
} // namespace abstract_squidget
