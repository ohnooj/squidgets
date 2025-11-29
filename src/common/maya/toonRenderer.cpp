#include "common/maya/toonRenderer.hpp"

namespace common {
namespace maya {
MObjectArray ToonRenderer::renderCurves(MStringArray outputGeom,
                                        bool worldSpace) {
  setupToonNode(outputGeom);
  MObjectArray ret = createToonCurves(worldSpace);

  MObject toonNode = toonNodePath.transform();
  MGlobal::deleteNode(toonNode);

  return ret;
}

void ToonRenderer::setupToonNode(MStringArray geomToRender) {
  MSelectionList candidateObjectList;
  for (MString s : geomToRender) {
    candidateObjectList.add(s);
  }
  MGlobal::setActiveSelectionList(candidateObjectList);

  // Create Toon Node.  Assume object to render is selected.
  MStatus status;
  status = MGlobal::executeCommand("AssignNewPfxToon");

  // Save toon node path.  Selection is on toon node
  MSelectionList sel_list;
  MGlobal::getActiveSelectionList(sel_list);
  sel_list.getDagPath(0, toonNodePath);

  //   The code below makes sure I render curves from the correct camera point.
  // For some reason (though calling AssignNewPfxToon in mel works), the
  // outlines render from some arbitrary point (i think the origin), so the
  // from these outlines are not from the correct view.
  //   I don't remember what the code below does. I think what I did was I moved
  // the camera to each object, rendered the curves, and then moved the camera
  // back.

  // Setup Toon Node's Camera Point for correct curve rendering.
  // Camera
  M3dView view = M3dView::active3dView();
  MDagPath camera;
  view.getCamera(camera);
  MFnTransform fnCamTransform(camera.transform());
  MTransformationMatrix camMatrix = fnCamTransform.transformation();
  MVector camTranslate = camMatrix.getTranslation(MSpace::kWorld);

  // Set Camera Point
  MPlugArray parr;
  MFnDagNode fnDagNode;
  fnDagNode.setObject(toonNodePath);
  fnDagNode.getConnections(parr);

  // Remove hierarchy shapes.
  for (MPlug p : parr) {
    if (p.isDestination()) {
      MFnDagNode fnPlugNode(p.source().node());
      MDagPath path;
      fnPlugNode.getPath(path);
      MFnTransform fnTrans(path.transform());

      if (geomToRender.indexOf(fnTrans.partialPathName()) == -1) {
        status = MGlobal::executeCommand(MString("disconnectAttr ") +
                                         p.source().name() + " " + p.name());
      }
    }
  }

  MObject cpAttr;
  MPlug cpPlug;
  cpAttr = fnDagNode.attribute("cameraPointX", &status);
  cpPlug = fnDagNode.findPlug(cpAttr, false);
  cpPlug.setDouble(camTranslate.x);

  cpAttr = fnDagNode.attribute("cameraPointY", &status);
  cpPlug = fnDagNode.findPlug(cpAttr, false);
  cpPlug.setDouble(camTranslate.y);

  cpAttr = fnDagNode.attribute("cameraPointZ", &status);
  cpPlug = fnDagNode.findPlug(cpAttr, false);
  cpPlug.setDouble(camTranslate.z);
  return;
}

MObjectArray ToonRenderer::createToonCurves(bool worldSpace) {
  /**
   * @brief Renders curves from outputGeom from toon node into 3D space.
   */
  MObjectArray ret;
  MStatus status;

  // Get render lines info.
  MFnPfxGeometry fnSet(toonNodePath);
  MRenderLineArray mainLines, leafLines, flowerLines;
  status =
      fnSet.getLineData(mainLines, leafLines, flowerLines, true, false, false,
                        false, false, false, false, false, worldSpace);

  // Create NURBS from curve geometry
  for (int i = 0; i < mainLines.length(); ++i) {
    MRenderLine renderLine = mainLines.renderLine(i, &status);
    MVectorArray renderPoints = renderLine.getLine();

    MPointArray curvePoints;
    for (unsigned j = 0; j < renderPoints.length(); ++j) {
      curvePoints.append(MPoint(renderPoints[j]));
    }

    int deg = 0;
    MObject curve = maya::CreateNurbsCurve(curvePoints, deg);
    MFnDagNode fnDagNode(curve);
    MDagPath dagPath;
    fnDagNode.getPath(dagPath);
    ret.append(dagPath.node());
  }

  return ret;
}
} // namespace maya
} // namespace common
