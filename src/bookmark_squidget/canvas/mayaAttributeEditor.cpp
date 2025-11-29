#include "bookmark_squidget/canvas/mayaAttributeEditor.hpp"

namespace common {
namespace maya {

MayaAttributeEditor::MayaAttributeEditor() {}

MayaAttributeEditor::~MayaAttributeEditor() {}

void MayaAttributeEditor::setObjectList(MSelectionList &objectList) {
  transform_paths = objectList;
}

AttrMap MayaAttributeEditor::getObjectAttributeMap() {
  return deformableObjectAttrMap;
}

void MayaAttributeEditor::createAttrMap() {
  MStatus status;

  // Find nodes with editable attributes: creationg and transform
  MObjectArray objectNodes;
  // status = getSourceCreationNodes(sourceCreationNodes);
  // status = getSourceTransformNodes(sourceTransformNodes);
  status = getObjectNodes(objectNodes);

  // Create map of objects : Attributes
  status = getObjectAttrMap(objectNodes, deformableObjectAttrMap);
}

MStatus MayaAttributeEditor::getObjectNodes(MObjectArray &objectNodes) {
  MDagPath path;
  for (unsigned i = 0; i < transform_paths.length(); ++i) {
    transform_paths.getDagPath(i, path);
    objectNodes.append(path.node());
  }

  return MStatus::kSuccess;
}

MStatus
MayaAttributeEditor::getObjectAttrMap(const MObjectArray &candidateObjects,
                                      AttrMap &sourceAttrMap) {
  /**
   * @brief This method searches for specific candidate attributes that can
   * deform our object.  To search for attributes, we search through the node
   * graph for the following nodes:
   * - Transform Nodes
   * - Object creation Nodes
   * For each object node
   *
   * Ofcourse there are many attributes to search for, but we restrict to
   * specific deformation attributes.
   */
  MStatus status;
  MDagPath path;
  MFnDagNode fnDagNode;

  for (MObject obj : candidateObjects) {
    fnDagNode.setObject(obj);
    fnDagNode.getPath(path);

    MapKey key = MString(path.fullPathName()).asChar();
    MPlugArray val = getEditableAttributes(obj);
    sourceAttrMap.emplace(key, val);
  }

  // Printing: plug name and vals
  for (auto &pair : sourceAttrMap) {
    auto key = pair.first;
    MPlugArray arr = sourceAttrMap.at(key);

    for (MPlug plug : arr) {
      float val;
      plug.getValue(val);
      // cout << "KeyVal: " << plug.name() << " : " << val << endl;
    }
  }

  return status;
}

MPlugArray MayaAttributeEditor::getEditableAttributes(MObject &obj) {
  MPlugArray retArr;
  MDagPath dagPath;
  MFnDagNode fnDagNode;

  std::unordered_set<std::string> validTransformAttrs;
  validTransformAttrs.emplace("translateX");
  validTransformAttrs.emplace("translateY");
  validTransformAttrs.emplace("translateZ");
  // validTransformAttrs.emplace("scaleX");
  // validTransformAttrs.emplace("scaleY");
  // validTransformAttrs.emplace("scaleZ");

  if (obj.hasFn(MFn::Type::kTransform)) {
    MFnTransform transform(obj);
    fnDagNode.setObject(obj);
    // MPlugArray connectionAttr;
    // transform.getConnections(connectionAttr);

    unsigned count = transform.attributeCount();
    for (unsigned i = 0; i < count; ++i) {
      MObject attr = transform.attribute(i);
      MFnAttribute fnAttr(attr);
      MPlug plug = fnDagNode.findPlug(attr, true);

      bool isValidPlug =
          fnAttr.isWritable() && fnAttr.isReadable() && !plug.isNull();
      bool isTransformAttr = validTransformAttrs.find(fnAttr.name().asChar()) !=
                             validTransformAttrs.end();
      if (isValidPlug && isTransformAttr) {
        retArr.append(plug);
      }
    }
  }

  return retArr;
}

void MayaAttributeEditor::updateAttrMap(AttrValueMap attr_val_map) {
  for (auto pair : attr_val_map) {
    auto key = pair.first;
    Eigen::VectorXd plugUpdatesValues = attr_val_map.at(key);
    MPlugArray arr = deformableObjectAttrMap.at(key);

    for (long i = 0; i < plugUpdatesValues.size(); ++i) {
      MPlug plug = arr[i];
      double update_val = plugUpdatesValues[i];
      plug.setValue(update_val);
    }
  }
}

AttrValueMap MayaAttributeEditor::getAttrValueMap() {
  AttrValueMap ret;
  for (auto pair : deformableObjectAttrMap) {
    auto key = pair.first;
    MPlugArray arr = deformableObjectAttrMap.at(key);

    Eigen::VectorXd val(arr.length());
    for (unsigned i = 0; i < arr.length(); ++i) {
      MPlug plug = arr[i];
      val[i] = plug.asDouble();
    }
    ret.emplace(key, val);
  }

  return ret;
}

// Unused methods
MStatus
MayaAttributeEditor::getSourceCreationNodes(MObjectArray &creationNodes) {
  /**
   * @brief Find construction nodes for shapes.
   *
   * @param shapeNodes
   */

  // cout << "getSourceCreationNodes" << endl;
  MDagPath dagPath;
  for (unsigned i = 0; i < transform_paths.length(); ++i) {
    transform_paths.getDagPath(i, dagPath);
    dagPath.extendToShapeDirectlyBelow(0);

    MFnDagNode fnDagNode(dagPath);
    MPlugArray plugArr;
    fnDagNode.getConnections(plugArr);

    for (MPlug plug : plugArr) {
      MFnAttribute fnAtt(plug.attribute());
      bool isInMesh = plug.asMObject().apiType() == MFn::Type::kMeshData;
      if (isInMesh) {
        MPlug source = plug.source();
        MObject sourceShape = source.node();
        creationNodes.append(sourceShape);
      }
    }
  }

  return MStatus::kSuccess;
}

MStatus
MayaAttributeEditor::getSourceTransformNodes(MObjectArray &transformNodes) {
  /**
   * @brief Get transform nodes from shape nodes.
   */

  MDagPath dagPath;
  for (unsigned i = 0; i < transform_paths.length(); ++i) {
    transform_paths.getDagPath(i, dagPath);
    transformNodes.append(dagPath.transform());
  }
  return MStatus::kSuccess;
}

MStatus MayaAttributeEditor::getToonRenderNode(MObject &toonNode) {
  /**
   * @brief Get toon render node from scene.
   */
  MStatus status;
  MSelectionList sList;
  MDagPath mDagPath;
  MObject mComponent;

  MGlobal::getActiveSelectionList(sList);
  MItSelectionList toonIter(sList, MFn::kPfxToon, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  status = toonIter.getDagPath(mDagPath, mComponent);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  toonNode = mDagPath.node(&status);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  return MStatus::kSuccess;
}

MStatus MayaAttributeEditor::getToonRenderSourceShapes(MObject &toonNode,
                                                       MObjectArray &retArray) {
  /**
   * @brief Get input shape nodes feeding into the toon render node.
   * From here we can find the transforms and creation attributes.
   *
   */
  MFnDagNode fnDagNode(toonNode);

  MPlugArray toonConnections;
  fnDagNode.getConnections(toonConnections);

  for (MPlug &p : toonConnections) {
    if (p.isDestination() && p.asMObject().apiType() == MFn::Type::kMeshData) {
      MPlug source = p.source();
      MObject sourceShape = source.node();
      retArray.append(sourceShape);
    }
  }

  return MStatus::kSuccess;
}
} // namespace maya
} // namespace common
