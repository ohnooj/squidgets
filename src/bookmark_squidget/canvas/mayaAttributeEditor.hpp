#ifndef MAYA_ATTRIBUTE_EDITOR_HPP
#define MAYA_ATTRIBUTE_EDITOR_HPP

#include <string>
#include <unordered_set>

#include <maya/MDagPath.h>
#include <maya/MFn.h>
#include <maya/MGlobal.h>
#include <maya/MObjectArray.h>
#include <maya/MObjectHandle.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MSelectionList.h>

#include <maya/MItSelectionList.h>

#include <maya/MFnAttribute.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnTransform.h>

#include "common/util/definitions.hpp"

namespace common {
namespace maya {

namespace {
using MapKey = util::MapKey;
using AttrMap = util::AttrMap;
using AttrValueMap = util::AttrValueMap;
} // namespace

/**
 * @brief Captures maya scene objects and adjusts attributes.
 *
 */
class MayaAttributeEditor {
public:
  MayaAttributeEditor();
  ~MayaAttributeEditor();

  void setObjectList(MSelectionList &objectList);
  void createAttrMap();

  void updateAttrMap(AttrValueMap attrValMap);
  AttrValueMap getAttrValueMap();

  AttrMap getObjectAttributeMap();

private:
  MStatus getToonRenderNode(MObject &toonNode);
  MStatus getToonRenderSourceShapes(MObject &toonNode, MObjectArray &retArray);

  MStatus getObjectNodes(MObjectArray &objectNodes);

  MStatus getSourceCreationNodes(MObjectArray &creationNodes);
  MStatus getSourceTransformNodes(MObjectArray &transformNodes);
  MStatus getObjectAttrMap(const MObjectArray &candidateObjects,
                           AttrMap &sourceAttrMap);
  MPlugArray getEditableAttributes(MObject &obj);

  MSelectionList transform_paths;
  AttrMap deformableObjectAttrMap;
};

} // namespace maya
} // namespace common

#endif
