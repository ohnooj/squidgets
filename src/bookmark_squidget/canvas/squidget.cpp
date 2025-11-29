#include "bookmark_squidget/canvas/squidget.hpp"

namespace bookmark_squidget {
namespace canvas {
Squidget::Squidget() {
  attributeConfigManager = std::make_shared<canvas::AttributeConfigManager>();
  elementContainer = std::make_shared<canvas::interface::ElementContainer>();
  t = -1;
  isSingleStroke = false;
}

Squidget::Squidget(MDagPath canvasPath) {
  attributeConfigManager = std::make_shared<canvas::AttributeConfigManager>();
  elementContainer =
      std::make_shared<canvas::interface::ElementContainer>(canvasPath);
  t = -1;
  isSingleStroke = false;
}

Squidget::~Squidget() {
  deactivate();
  selectionList.clear();
}

void Squidget::activate() {
  if (isSingleStroke) {

  } else {
    setDrivenConfigs();
  }
}

void Squidget::deactivate() {
  if (isSingleStroke) {

  } else {
    clearDrivenKeyframes();
  }
}

// =============================================================================
// CREATE MODE
// =============================================================================
void Squidget::addStepElement(
    std::shared_ptr<interface::StepElement> stepElement) {
  elementContainer->AddStepElement(stepElement);
}

void Squidget::addAttributeConfig(
    std::shared_ptr<AttrConfigMap> attrConfigMap) {
  attributeConfigManager->addConfig(attrConfigMap);
}

// =============================================================================
// EDIT MODE SETUP
// =============================================================================
void Squidget::setDrivenConfigs() {
  cout << "Squidget::setDrivenConfigs()" << endl;
  /**
   * @brief Called when edit mode is toggled.  Will call setDrivenKeys
   *
   * Get stroke steps and attribute configurations and sets driven keys with
   * stroke path value and driven values.
   */
  MStatus s;
  auto steps = elementContainer->m_stepElementKey;
  auto configs = attributeConfigManager->getChangedConfigs(isSingleStroke);

  MFnDagNode fnDagNode(elementContainer->pathShape);
  MPlug squidget_t_plug =
      fnDagNode.findPlug(elementContainer->squidget_t, true);

  for (size_t i = 0; i < steps.size(); ++i) {
    double driver_val = steps[i];
    AttributeConfigStep acs = configs[i];
    // if (isSingleStroke) {
    //   acs = configs[1];
    // }
    std::shared_ptr<AttrConfigMap> acm = acs.attrValueMap;
    squidget_t_plug.setDouble(driver_val);

    for (auto pair : *acm) {
      MString k(pair.first.c_str());
      double val = pair.second;

      std::ostringstream ss;
      ss << "cmds.setAttr('" << k << "', " << pair.second << ")" << endl;
      ;
      std::string cmd = ss.str();
      MGlobal::executePythonCommand(cmd.c_str());

      MDagPath pathShape = elementContainer->pathShape;
      MString pathDriver = pathShape.fullPathName() + "." + "squidget_t";

      // Example: setDrivenKeyframe -currentDriver PathShape.squidget_t
      // L_leg_IK_ctrl.translateZ;
      ss.str("");
      ss.clear();
      ss << "setDrivenKeyframe";
      ss << " -itt \"linear\" -ott \"linear\"";
      ss << " -currentDriver " << pathDriver;
      ss << " " << k << ";";
      cmd = ss.str();
      s = MGlobal::executeCommand(MString(cmd.c_str()));

      drivenKeyframes.append(k);
    }
  }
  setTValue(0);
}

void Squidget::clearDrivenKeyframes() {
  for (MString kf : drivenKeyframes) {
    std::ostringstream ss;
    ss << "cmds.cutKey(\"" << kf << "\")";
    MGlobal::executePythonCommand(ss.str().c_str());
  }

  drivenKeyframes.clear();
}

void Squidget::setDiscrete() {
  auto configs = attributeConfigManager->getChangedConfigs(isSingleStroke);

  if (t == 1) {
    AttributeConfigStep acs = configs[1];
    std::shared_ptr<AttrConfigMap> acm = acs.attrValueMap;

    for (auto pair : *acm) {
      MString k(pair.first.c_str());
      double val = pair.second;

      std::ostringstream ss;
      ss << "cmds.setAttr('" << k << "', " << pair.second << ")" << endl;
      ;
      std::string cmd = ss.str();
      MGlobal::executePythonCommand(cmd.c_str());
    }
  }
}

// =============================================================================
// EDIT MODE POINT
// =============================================================================
void Squidget::ReadPoint(MPoint localPoint) {
  double t = elementContainer->UpdateValueAlongPath(localPoint);
  lastStroke = elementContainer->getGhostShape();
}

bool Squidget::setTValue(double t) {
  elementContainer->setTValue(t);
  this->t = t;
  return true;
}

// =============================================================================
// EDIT MODE STROKE
// =============================================================================
void Squidget::SetValueFromStroke(const MPointArray localStroke,
                                  MPoint &pathPoint) {
  /**
   * @brief Creates stroke steps for the squidget and ties them to the
   * attribute configuration.
   *
   * How do I select the correct attribute configuration for a step element?
   */
  if (localStroke.length() == 0) {
    return;
  }

  // if stroke crosses path, use path/stroke intersection point as stroke
  bool crosses = checkStrokeCrossing(localStroke);
  double dist;
  if (crosses) {
    t = elementContainer->CheckPathCrossing(localStroke);
    dist = 0;
  } else {
    t = estimateStrokeValue(localStroke, dist);
  }

  cout << "t val: " << t << " dist: " << dist << endl;
  setTValue(t);
  pathPoint = elementContainer->getPointAlongPath(
      t * (elementContainer->m_stepElements.size() -
           1)); // parameter is not from 0 to 1...
  lastStroke = localStroke;

  if (isSingleStroke) {
    setTValue(1);
    setDiscrete();
  }
}

double Squidget::estimateStrokeValue(const MPointArray localPts, double &dist) {
  MPointArray reverseLocalPoints;
  for (int i = localPts.length() - 1; i >= 0; --i)
    reverseLocalPoints.append(localPts[i]);

  double revDist;
  double t = elementContainer->queryStrokeValue(localPts, dist);
  return t;
}

void Squidget::setSingleStrokeAttributes() {
  attributeConfigManager->interpolate(1);
}

bool Squidget::checkStrokeCrossing(MPointArray localPts) {
  return elementContainer->CheckPathCrossing(localPts) != -1;
}

MObject Squidget::getPathShape() { return elementContainer->pathShape.node(); }

} // namespace canvas
} // namespace bookmark_squidget