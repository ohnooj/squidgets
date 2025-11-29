#include "bookmark_squidget/canvas/squidgetCanvas.hpp"

namespace bookmark_squidget {
namespace canvas {
SquidgetCanvas::SquidgetCanvas() {}

SquidgetCanvas::SquidgetCanvas(MDagPath canvasPath) {
  _canvasPath = canvasPath;
  id = common::util::generateRandomID(2);

  // Relies on selected objects
  rest_attrConfigMap = std::make_shared<AttrConfigMap>();
  setRestAttributeConfig();
}

SquidgetCanvas::~SquidgetCanvas() {
  MObject pt = _canvasPath.transform();
  MGlobal::deleteNode(pt);
}

void SquidgetCanvas::activate() {
  // Find indices not in squidgets
  MIntArray indices;
  for (size_t i = 0; i < m_stepElements.size(); ++i) {
    bool isUsed = false;

    auto stepElement = m_stepElements[i];
    for (auto &s : m_squidgets) {
      auto it =
          std::find(s->elementContainer->m_stepElements.begin(),
                    s->elementContainer->m_stepElements.end(), stepElement);
      if (it != s->elementContainer->m_stepElements.end()) {
        isUsed = true;
        break;
      }
    }

    if (!isUsed) {
      indices.append(i);
    }
  }

  for (int i : indices) {
    MIntArray arr;
    arr.append(i);
    addSquidget(arr);
  }
}

void SquidgetCanvas::deactivate() {
  size_t i = 0;
  while (i < m_squidgets.size()) {
    bool isSingleStroke = m_squidgets[i]->getSingleStroke();

    if (isSingleStroke) {
      m_squidgets.erase(m_squidgets.begin() + i);
    } else {
      i++;
    }
  }
}

int SquidgetCanvas::restAttrSize() { return rest_attrConfigMap->size(); }

// =============================================================================
//                            Get Attribute Config
// =============================================================================
void SquidgetCanvas::setRestAttributeConfig() {
  MGlobal::getActiveSelectionList(restSelectionList);

  MStringArray goodPlugs;
  getSelectedAttributePlugs(goodPlugs);
  if (goodPlugs.length() == 0) {
    MGlobal::displayError("No animatable attributes selected.");
    return;
  }

  rest_attrConfigMap = getAttributeConfig(goodPlugs);

  MStringArray names;
  restSelectionList.getSelectionStrings(names);
}

bool SquidgetCanvas::isValid() {
  bool canvasPathValid = _canvasPath.isValid();
  size_t i = 0;
  while (i < m_stepElements.size()) {
    if (m_stepElements[i]->isValid()) {
      i++;
    } else {
      m_stepElements.erase(m_stepElements.begin() + i);
      m_attributeConfigMaps.erase(m_attributeConfigMaps.begin() + i);
    }
  }

  return _canvasPath.isValid();
}

void SquidgetCanvas::setVisibility(bool isVisible) {
  MFnDagNode fnDagNode(_canvasPath.transform());
  int isVisibleInt = isVisible ? 1 : 0;
  std::ostringstream ss;
  ss << "cmds.setAttr('" << fnDagNode.fullPathName() << ".visibility', ";
  ss << isVisibleInt << ")";
  MGlobal::executePythonCommand(ss.str().c_str());
}

void SquidgetCanvas::highlight(bool highlight) {
  std::ostringstream ss;
  if (highlight) {
    ss << "cmds.sets(\"" << _canvasPath.fullPathName()
       << "\", e=True, forceElement=\"canvasHighlightColorSG\");";
    MGlobal::executePythonCommand(ss.str().c_str());
  } else {
    ss << "cmds.sets(\"" << _canvasPath.fullPathName()
       << "\", e=True, forceElement=\"canvasNormalColorSG\");";
    MGlobal::executePythonCommand(ss.str().c_str());
  }
}

bool SquidgetCanvas::isVisible() {
  MFnDagNode fnDagNode(_canvasPath.transform());
  std::ostringstream ss;
  ss << "cmds.getAttr('" << fnDagNode.fullPathName() << ".visibility')";
  MString res = MGlobal::executePythonCommandStringResult(ss.str().c_str());
  return (res == "True");
}

// =============================================================================
//                            CreateMode Stroke
// =============================================================================
void SquidgetCanvas::HandleCreateModeStroke(MPointArray screenStroke,
                                            MSelectionList selList,
                                            bool &shouldDelete) {
  selectionList = selList;

  bool isStrokeLong = screenStroke.length() >= 5;
  if (!isStrokeLong)
    return;

  bool isVis = isVisible();

  MPointArray worldStroke;
  bool intersectsCanvas = getIntersectingWorldPoints(screenStroke, worldStroke);
  if (!intersectsCanvas)
    return;
  MPointArray localStroke = localizePointsToCanvas(worldStroke);

  MIntArray crossingElementIndices = checkElementCrossing(localStroke);
  bool isCrossingElements = crossingElementIndices.length() > 0;
  MIntArray crossingPathIndicies = checkPathCrossing(localStroke);
  bool isCrossingSquidgetPath = crossingPathIndicies.length() > 0;
  bool isSelfXing = checkSelfXing(localStroke);

  if (isSelfXing) {
    bool hasSquidget = m_squidgets.size() > 0;
    bool hasElements = m_stepElements.size() > 0;
    shouldDelete = !(hasSquidget || hasElements);
    return;
  }

  if (isCrossingSquidgetPath) {
    // Crosses only Squidget Path
    if (crossingPathIndicies.length() == 1) {
      return;
    }

    // Must cross twice to delete
    if (crossingPathIndicies[0] == crossingPathIndicies[1]) {
      std::string first = m_squidgets[crossingPathIndicies[0]]
                              ->elementContainer->m_stepElements[0]
                              ->id;
      std::string last =
          m_squidgets[crossingPathIndicies[0]]
              ->elementContainer
              ->m_stepElements[m_squidgets[crossingPathIndicies[0]]
                                   ->elementContainer->m_stepElements.size() -
                               1]
              ->id;

      m_squidgets.erase(m_squidgets.begin() + crossingPathIndicies[0]);
    }
  } else if (isCrossingElements && !isCrossingSquidgetPath) {
    // Crosses only Elements
    if (crossingElementIndices.length() == 1) {
      return;
    }

    // Must scratch only one element.
    int toRemove = -1;
    bool crossesTwice = crossingElementIndices[0] == crossingElementIndices[1];
    bool noConnected =
        m_stepElements[crossingElementIndices[0]].use_count() == 1;
    toRemove = crossingElementIndices[0];

    if (crossesTwice &&
        noConnected) { // Make sure element isn't part of squidget
      removeStepElement(toRemove);
    } else {
      bool allElementsNaked = true;
      for (int i : crossingElementIndices) {
        if (m_stepElements[i].use_count() > 1) {
          allElementsNaked = false;
        }
      }
      if (allElementsNaked) {
        addSquidget(crossingElementIndices);
      }
    }
  } else if (!isCrossingElements && !isCrossingSquidgetPath) {
    addStepElement(localStroke);
    addAttributeConfig();
  } else {
    MGlobal::displayWarning("Double cross over a line to delete it.");
  }
}

void SquidgetCanvas::addStepElement(MPointArray localStroke) {
  m_stepElements.emplace_back(createStepElement(localStroke));
}

void SquidgetCanvas::removeStepElement(int index) {
  m_stepElements.erase(m_stepElements.begin() + index);
  m_attributeConfigMaps.erase(m_attributeConfigMaps.begin() + index);
}

std::shared_ptr<AttrConfigMap>
SquidgetCanvas::getAttributeConfig(MStringArray plugs) {
  MStatus res;
  AttrConfigMap acm;
  for (MString s : plugs) {
    std::ostringstream ss;
    ss << "cmds.getAttr('" << s.asChar() << "')";
    std::string cmd = ss.str();
    MString val = MGlobal::executePythonCommandStringResult(cmd.c_str(), false,
                                                            false, &res);
    if (res != MStatus::kSuccess) {
      continue;
    }

    double d = val.asDouble();
    if (val == "True") {
      d = 1;
    } else if (val == "false") {
      d = 0;
    }

    acm.emplace(s.asChar(), d);
  }

  return std::make_shared<AttrConfigMap>(acm);
}

void SquidgetCanvas::addAttributeConfig() {
  MStringArray plugNames;
  getSelectedAttributePlugs(plugNames);

  auto new_acm = getAttributeConfig(plugNames);

  // Calculate intersection w/ rest config saved into a new map
  AttrConfigMap intersect;
  for (auto it = new_acm->begin(); it != new_acm->end(); ++it) {
    auto restIt = rest_attrConfigMap->find(it->first);
    if (restIt != rest_attrConfigMap->end()) {
      if (it->second != restIt->second) {
        intersect.emplace(it->first, it->second);
      }
    }
  }
  m_attributeConfigMaps.push_back(new_acm);
}

void SquidgetCanvas::addSquidget(MIntArray indices) {
  std::shared_ptr<Squidget> sq = std::make_shared<Squidget>(_canvasPath);

  if (indices.length() == 1) { // Discrerte Squidget
    int i = indices[0];
    sq->addStepElement(m_stepElements[i]);
    sq->addAttributeConfig(rest_attrConfigMap);
    sq->addStepElement(m_stepElements[i]);
    sq->addAttributeConfig(m_attributeConfigMaps[i]);
    sq->setSingleStroke(true);
  } else {
    std::string first_id = m_stepElements[indices[0]]->id;
    std::string last_id = m_stepElements[indices[indices.length() - 1]]->id;
    for (int i : indices) {
      sq->addStepElement(m_stepElements[i]);
      sq->addAttributeConfig(m_attributeConfigMaps[i]);
    }
  }
  sq->activate();
  m_squidgets.emplace_back(sq);
}

// =============================================================================
//                            EditMode Stroke
// =============================================================================
bool SquidgetCanvas::HandleEditModeStroke(const MPointArray screenStroke,
                                          MPoint &pathPoint) {
  // Query the correct Squidget.
  bool isStrokeLong = screenStroke.length() >= 5;
  if (!isStrokeLong) {
    return false;
  }

  MPointArray worldStroke;
  bool intersectsCanvas = getIntersectingWorldPoints(screenStroke, worldStroke);
  if (!intersectsCanvas) {
    return false;
  }

  MPointArray localStroke = localizePointsToCanvas(worldStroke);
  _activeSquidgetIndex = activeSquidgetIndex(localStroke);
  if (_activeSquidgetIndex < 0) {
    return false;
  }

  m_squidgets[_activeSquidgetIndex]->SetValueFromStroke(localStroke, pathPoint);
  return true;
}

void SquidgetCanvas::ReadPoint(MPoint screenPoint) {
  /**
   * Prior: activeSquidgetIndex is set.
   */
  MPointArray screenPoints;
  screenPoints.append(screenPoint);

  MPointArray worldPoints;
  // Comment out for now
  bool intersectsElement =
      getIntersectingWorldPoints(screenPoints, worldPoints);
  if (!intersectsElement)
    return;

  // index is valid because of callback
  MPoint localPt = localizePointsToCanvas(worldPoints)[0];
  if (_activeSquidgetIndex < 0 ||
      _activeSquidgetIndex >= (int)m_squidgets.size()) {
    // For some reason i get an index out of bounds sometimes.
    unhandleSquidget();
  } else {
    m_squidgets[_activeSquidgetIndex]->ReadPoint(localPt);
  }
}

// =============================================================================
//                            Extra
// =============================================================================
void SquidgetCanvas::getSelectedAttributePlugs(MStringArray &arr) {
  // cout << "getSelectedAttributePlugs()" << endl;
  MGlobal::setActiveSelectionList(restSelectionList);

  std::ostringstream ss;
  ss << "def getGoodPlugs():" << endl;
  ss << "    sss = cmds.ls(sl=True, l=True)" << endl;
  ss << "    " << endl;
  ss << "    relatives = []" << endl;
  ss << "    while len(sss) > 0:" << endl;
  ss << "        curr = sss.pop()" << endl;
  ss << "        relatives.append(curr)" << endl;
  ss << "        if cmds.getAttr(f'{curr}.visibility'):" << endl;
  ss << "            curr_rels = cmds.listRelatives(curr, f=True, "
        "typ='transform') or []"
     << endl;
  ss << "            sss.extend(curr_rels)" << endl;

  // ss << "    sss = cmds.ls(sl=True, l=True)"<<endl;
  // ss << "    relatives = cmds.listRelatives(sss, ad=True, f=True,
  // typ='transform') or []"<<endl;
  ss << "    relatives.extend(sss)" << endl;
  ss << "    relatives = [r for r in relatives if 'Constraint' not in "
        "cmds.nodeType(r)]"
     << endl;
  ss << "    relatives = [r for r in relatives if "
        "cmds.getAttr(f'{r}.visibility')]"
     << endl;
  ss << "    animatableAttrs = []" << endl;
  ss << "    [animatableAttrs.extend(cmds.listAnimatable(v) or []) for v in "
        "relatives]"
     << endl;

  ss << "    user_defined = []" << endl;
  ss << "    if len(relatives) > 0:" << endl;
  ss << "        user_defined=[]" << endl;
  ss << "        for r in relatives:" << endl;
  ss << "            user_attrs = cmds.listAttr(r, userDefined=True, "
        "scalar=True) or []"
     << endl;
  ss << "            user_attrs_long = [f'{r}.{att}' for att in user_attrs]"
     << endl;
  ss << "            user_defined.extend(user_attrs_long)" << endl;
  ss << "        history = cmds.listHistory(relatives)" << endl;
  ss << "        bs=cmds.ls(history, type='blendShape')" << endl;
  ss << "        bs_anim=cmds.listAnimatable(bs) or []" << endl;

  ss << "    animatableAttrs.extend(user_defined)" << endl;
  ss << "    animatableAttrs.extend(bs_anim)" << endl;

  ss << "    retAttrs = []" << endl;
  ss << "    for att in animatableAttrs:" << endl;
  ss << "        isSettableKeyableUnlocked = cmds.getAttr(att, se=True, "
        "l=False, k=True)"
     << endl;
  ss << "        isDest = cmds.connectionInfo(att, id=True)" << endl;
  ss << "        isExactDest = cmds.connectionInfo(att, ied=True)" << endl;
  ss << "        isLocked = cmds.connectionInfo(att, il=True)" << endl;
  ss << "        isJoint = cmds.nodeType(att) == 'joint'" << endl;
  ss << "        isGood = isJoint or isSettableKeyableUnlocked and (not "
        "(isDest or isExactDest or isLocked))"
     << endl;
  ss << "        if (isGood):" << endl;
  ss << "            retAttrs.append(att)" << endl;
  ss << "    return retAttrs" << endl;

  MString cmd = ss.str().c_str();
  MGlobal::executePythonCommand(cmd);

  MString goodPlugs =
      MGlobal::executePythonCommandStringResult("getGoodPlugs()");
  goodPlugs.substitute("'", "");
  goodPlugs.substitute(" ", "");
  goodPlugs = goodPlugs.substring(1, goodPlugs.length() - 2);

  cout << "Good plugs: " << goodPlugs << endl;
  arr.clear();
  goodPlugs.split(',', arr);
}

bool SquidgetCanvas::getIntersectingWorldPoints(const MPointArray screenPts,
                                                MPointArray &worldPts) {
  bool allIntersects = true;
  for (MPoint pt : screenPts) {
    MPoint worldPt;
    bool intersects = getIntersectingWorldPoint(pt, worldPt);

    worldPts.append(worldPt);
    allIntersects = allIntersects && intersects;
  }

  return allIntersects;
}

bool SquidgetCanvas::getIntersectingWorldPoint(const MPoint screenPt,
                                               MPoint &worldPt) {
  MFnNurbsSurface surfaceFn;
  surfaceFn.setObject(_canvasPath);

  MPoint worldCamOrigin;
  MVector worldCamVec;
  M3dView view = M3dView::active3dView();
  view.viewToWorld(screenPt.x, screenPt.y, worldCamOrigin, worldCamVec);
  // cout << worldCamOrigin << " " << worldCamVec << endl;
  double u, v;
  return surfaceFn.intersect(worldCamOrigin, worldCamVec, u, v, worldPt,
                             kMFnNurbsEpsilon, MSpace::kWorld);
}

int SquidgetCanvas::activeSquidgetIndex(MPointArray localPts) {
  // cout << "SquidgetCanvas::activeSquidgetIndex()" << endl;
  if (_activeSquidgetIndex < 0) {
    int index = _getClosestSquidgetIndex(localPts);
    if (index < 0) {
      // cout << "No Squidgets" << endl;
      return -1;
    }
    // cout << "Squidget Index: " << index << endl;
    _activeSquidgetIndex = index;
  }
  return _activeSquidgetIndex;
}

int SquidgetCanvas::_getClosestSquidgetIndex(MPointArray localPts) {
  // Only get a squidget if the stroke crosses the squidget path.
  double minDist = std::numeric_limits<double>::max();
  int minIndex = -1;

  for (size_t i = 0; i < m_squidgets.size(); ++i) {
    // Should immediately use squidget with cross over.
    bool crosses = m_squidgets[i]->checkStrokeCrossing(localPts);
    if (crosses) {
      return i;
    }

    double dist;
    double t_val = m_squidgets[i]->estimateStrokeValue(localPts, dist);
    // if (t_val > 0 && t_val < 1. && dist < minDist) {
    if (dist < minDist) {
      minDist = dist;
      minIndex = i;
    }
  }
  return minIndex;
}

// =============================================================================
//                            Checks
// =============================================================================
MIntArray SquidgetCanvas::checkPathCrossing(MPointArray localStroke) {
  // cout << "SquidgetCanvas::CheckPathCrossing()" << endl;
  MPointArray stroke =
      common::util::swapPointArrayCoordinates(localStroke, 1, 2);
  MIntArray ret;
  for (size_t i = 0; i < m_squidgets.size(); i++) {
    auto squidget = m_squidgets[i];
    MFnNurbsCurve curveFn(squidget->getPathShape());
    MPointArray cvs;
    curveFn.getCVs(cvs);
    MPointArray pathShape = common::util::swapPointArrayCoordinates(cvs, 1, 2);

    int num_intersections = common::util::polylinesIntersect(stroke, pathShape);
    for (int j = 0; j < num_intersections; ++j) {
      ret.append(i);
    }
  }
  return ret;
}

MIntArray SquidgetCanvas::checkElementCrossing(MPointArray localStroke) {
  // cout << "SquidgetCanvas::checkElementCrossing()" << endl;
  MPointArray stroke =
      common::util::swapPointArrayCoordinates(localStroke, 1, 2);
  MIntArray ret;

  for (unsigned int i = 0; i < stroke.length() - 1; ++i) {
    MPointArray strokeSegment;
    strokeSegment.append(stroke[i]);
    strokeSegment.append(stroke[i + 1]);

    for (size_t j = 0; j < m_stepElements.size(); j++) {
      auto stepElement = m_stepElements[j];
      MPointArray stepShape =
          common::util::swapPointArrayCoordinates(stepElement->shape, 1, 2);

      int num_intersections =
          common::util::polylinesIntersect(strokeSegment, stepShape);
      for (int k = 0; k < num_intersections; ++k) {
        ret.append(j);
      }
    }
  }

  for (int i : ret) {
    MObject o = m_stepElements[i]->curve;
    MDagPath path;
    MFnDagNode fnDagNode(o);
    fnDagNode.getPath(path);
  }

  return ret;
}

bool SquidgetCanvas::checkSelfXing(MPointArray localStroke) {
  // cout << "SquidgetCanvas::checkSelfXing()" << endl;
  MPointArray stroke =
      common::util::swapPointArrayCoordinates(localStroke, 1, 2);

  for (unsigned int i = 0; i < stroke.length() - 1; ++i) {
    MPointArray a;
    a.append(stroke[i]);
    a.append(stroke[i + 1]);

    for (unsigned int j = i + 2; j < stroke.length() - 1; ++j) {
      MPointArray b;
      b.append(stroke[j]);
      b.append(stroke[j + 1]);

      int num_intersections = common::util::polylinesIntersect(a, b);
      if (num_intersections > 0) {
        return true;
      }
    }
  }

  return false;
}

void SquidgetCanvas::unhandleSquidget() { _activeSquidgetIndex = -1; }

MPointArray SquidgetCanvas::localizePointsToCanvas(MPointArray pts) {
  // This doesn't handle for non-global parent hierarchy.  Using cmds.xform

  MFnTransform fnTransform(_canvasPath.transform());

  std::ostringstream ss;
  ss << "cmds.xform('";
  ss << fnTransform.fullPathName();
  ss << "', q=True, m=True, ws=True)";
  MString m = MGlobal::executePythonCommandStringResult(ss.str().c_str());
  m = m.substring(1, m.length() - 1);

  MStringArray strVals;
  m.split(',', strVals);
  MDoubleArray vals;
  for (MString s : strVals) {
    vals.append(s.asDouble());
  }

  double dArr[4][4];
  int index = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      dArr[i][j] = vals[index++];
    }
  }
  MMatrix worldXForm(dArr);
  MMatrix canvasMatrix = worldXForm.inverse();

  MPointArray localStroke; // apply inverse transformation to points.
  for (MPoint p : pts) {
    localStroke.append(p * canvasMatrix);
  }
  return localStroke;
}

std::shared_ptr<interface::StepElement>
SquidgetCanvas::createStepElement(MPointArray pts) {
  std::shared_ptr<canvas::interface::StepElement> stepElement =
      std::make_shared<canvas::interface::StepElement>();

  MObject canvasObj = _canvasPath.transform();
  MObject curve = common::maya::CreateNurbsCurve(pts, 0, canvasObj);
  curve = common::util::resampleCurve(curve, 30);

  MFnDagNode fnDagNode;
  fnDagNode.setObject(curve);
  std::ostringstream ss;
  ss << "rename " << fnDagNode.fullPathName() << " ";
  ss << "StepElement";
  // cout << "Command: " << ss.str() << endl;
  MGlobal::executeCommand(ss.str().c_str());

  stepElement->shape = pts;
  stepElement->curve = curve;
  return stepElement;
}
} // namespace canvas
} // namespace bookmark_squidget
