#ifndef ABSTRACT_SQUIDGET_TOOL_HPP
#define ABSTRACT_SQUIDGET_TOOL_HPP

#include <Eigen/Dense>
#include <vector>

#include <maya/MDagPath.h>
#include <maya/MMatrix.h>
#include <maya/MPointArray.h>

#include <maya/MPxToolCommand.h>

#include "abstract_squidget/types.hpp"
#include "common/logging/logger.hpp"
#include "common/util/definitions.hpp"

namespace abstract_squidget {
namespace apps {
#define kCurveFlag "-c"
#define kCurveFlagLong "-curve"
#define kCandidateObjectsFlag "-o"
#define kCandidateObjectsFlagLong "-objects"
#define kPreSelectFlag "-ps"
#define kPreSelectFlagLong "-preselect"
#define kDeformFlag "-d"
#define kDeformFlagLong "-deform"
#define kTranslateOnlyFlag "-t"
#define kTranslateOnlyFlagLong "-translate"
#define kShadowCurveFlag "-sh"
#define kShadowCurveFlagLong "-shadowcurve"
#define kIsRigCurvesFlag "-r"
#define kIsRigCurvesFlagLong "-rigcurves"

namespace {
using AttrValueMap = common::util::AttrValueMap;
using PenPoints = common::util::PenPoints;

} // namespace
/**
 * AST: Abstract Squidget Tool
 *  - curve: [points] 2D screen space input stroke
 *  - objs: list[str] 3D candidate objects to choose frome
 *  - select: bool    whether to select objects.  Will return (obj, screen_space
 * silhouette)
 *  - deform: bool    whether to deform objects.
 *
 *
 * CurveTool does the following:
 * - Deforming object from attribute
 * - From a given stroke, look through all attributes that would deform
 * - Selects object attributes from stroke points
 * - deforms attributes.
 *
 * Given a set of deformable objects Mi with attributes Ai and a stroke S_c, and
 * generating toon curve T_j from M_i*A_i, find an attribute config A'_i that
 * generates toon curves T'_j that minimizes curve_dist(S_c, t'_j). In other
 * words, find an attribute matrix that renders toon curves that's closest to
 * our stroke. Have a regularization of A'_i - A_i (possibly a L1
 * regularization).
 *
 * Things to address:
 * - Create toon curve renders T_j for MObject M_i
 * - Compare curve distance of T_j to S_c.
 *      - Compare in viewport domain.
 * - Find creation and transform attributes.
 *      - Question: how do we recognize creation vs transform attributes?  Maybe
 * there's a hierarchy.
 * - Optimize iteratively to find deformation.
 */

class AbstractSquidgetTool : public MPxToolCommand {
public:
  AbstractSquidgetTool();
  ~AbstractSquidgetTool() override;
  static void *creator();
  static MSyntax newSyntax();

  // MPxToolCommand Methods ========================
  MStatus doIt(const MArgList &args) override;
  MStatus redoIt() override;
  MStatus undoIt() override;
  MStatus parseArgs(const MArgList &args);

  bool isUndoable() const override;
  bool hasSyntax() const override;
  // MSyntax syntax() const override;
  // bool isHistoryOn() const override;
  // MString commandString() const override;
  // MStatus setHistoryOn(bool historyOn) override;
  // MStatus setCommandString(const MString &commandString) override;
  MArgList getFinalCommand();
  MString getCallableCommandStr();
  EditFlags getEditFlags();
  MStatus cancel() override;
  MStatus finalize() override;

  // Getters and Setters ========================
  void setCurvePoints2D(MString curvePoints2D);
  void setCandidateObjs(MString candidateObjs);
  void setPreSelect(bool preSelect);
  void setDeform(bool deform);
  void setTranslateOnly(bool translateOnly);
  void setCurveToMatchStr(MString shadowCurve);
  void setIsRigCurves(bool rigCurves);

  MString getCurvePoints2D() { return curvePoints2DStr; }
  MString getCandidateObjs() { return candidateObjsStr; }
  bool getPreSelect() { return preSelect; }
  bool getDeform() { return deform; }
  bool getTranslateOnly() { return translateOnly; }
  MString getCurveToMatchStr() { return inCurveToMatchStr; }
  bool getIsRigCurves() { return isRigCurves; }

  void setLogger(std::shared_ptr<common::logging::Logger> logger);

  // Back to Context Methods ========================
  // void exportOutShadowCurve(MString &shadowCurve);
  // void exportContextShadowCurveWorld();
  // void exportContextEditedObject();

  void fromContextShadowCurveWorld(std::shared_ptr<MPointArray> s);
  void fromContextShadowCurveLocalStr(std::shared_ptr<MString> s);
  void fromContextEditedObject(std::shared_ptr<MString> editedObject);

  void exportEditedObjectToContext(MString obj);
  void exportShadowCurveLocalToContext(MPointArray arr);
  void exportShadowCurveWorldToContext(MPointArray arr, RenderMatrices rms);

private:
  MString buildOutShadowCurveStr(MPointArray osc);
  MPointArray buildOutShadowCurveWorld(MPointArray outShadowCurve,
                                       RenderMatrices renderMats);

  void printCommand(MArgList command);
  void saveForUndo();

  // Input Variables ========================
  MString curvePoints2DStr = "";
  MPointArray curvePoints2D;
  MString candidateObjsStr = "";
  MStringArray candidateObjs;
  MString inCurveToMatchStr = "";
  MPointArray inCurveToMatch;

  // Output Variables ========================
  MString editedObject = "";
  MString outShadowCurveStr = "";
  MPointArray outShadowCurveLocal;
  MPointArray outShadowCurveWorld;

  bool preSelect = false; // Shift
  bool deform = false;    // Ctrl
  bool translateOnly = false;
  bool isRigCurves = false;

  // Undo saved variables
  MTransformationMatrix oldMatrix;
  MPointArray oldVertices;
  MPointArray oldShadowCurveWorld;
  MString oldShadowCurveLocalStr;
  MString oldEditedObj;

  // IMPORTANT: This is the main selection group used for the tool.
  int mode = 1; // study_0
                // 1; // study_1
                // 2; // study_2

  // From Context
  std::shared_ptr<common::logging::Logger> logger;
  std::shared_ptr<MPointArray> shadowCurveWorld;
  std::shared_ptr<MString> shadowCurveLocalStr;
  std::shared_ptr<MString> editedObj;
};

} // namespace apps
} // namespace abstract_squidget

#endif