#include "abstract_squidget/deformer/renderMatrices.hpp"

#include <maya/M3dView.h>
#include <maya/MDagPath.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnTransform.h>
#include <maya/MQuaternion.h>
#include <maya/MVector.h>

#include "common/util/maya.hpp"

namespace abstract_squidget {
namespace deformer {

RenderMatrices getRenderMatrices(MObject geomObj) {
  /**
   * @brief Returns all the matrices within the render pipeline.
   * Object -> Local -> World -> Camera -> Projection -> Screen
   */
  MFnDagNode fnDagNode(geomObj);
  MDagPath geomPath;
  fnDagNode.getPath(geomPath);
  MFnTransform fnTransform(geomPath);

  M3dView view = M3dView::active3dView();
  int width = view.portWidth();
  int height = view.portHeight();
  MMatrix cameraMatrix, projMatrix;
  view.modelViewMatrix(cameraMatrix);
  view.projectionMatrix(projMatrix);
  double screenArray[4][4] = {
      {((double)width) / 2, 0, 0, ((double)width) / 2},
      {0, ((double)height) / 2, 0, ((double)height) / 2},
      {0, 0, 1, 0},
      {0, 0, 0, 1}};
  cameraMatrix = cameraMatrix.transpose();
  projMatrix = projMatrix.transpose(); // Cause maya representation
  MMatrix screenMatrix(screenArray);
  // MMatrix worldMatrix = fnTransform.transformationMatrix().transpose();

  MMatrix localMatrix = fnTransform.transformation().asMatrix().transpose();

  MTransformationMatrix worldMat = fnTransform.transformation();
  MVector worldTrans;
  MQuaternion worldRot;
  worldTrans = fnTransform.getTranslation(MSpace::kWorld);
  fnTransform.getRotationQuaternion(worldRot.x, worldRot.y, worldRot.z,
                                    worldRot.w, MSpace::kWorld);
  worldMat.setRotationQuaternion(worldRot.x, worldRot.y, worldRot.z,
                                 worldRot.w);
  worldMat.setTranslation(worldTrans, MSpace::kWorld);

  // World = Parent * Local
  // World * Local^-1 = Parent
  // This is actually parent matrix
  MMatrix worldMatrix = worldMat.asMatrix().transpose();
  worldMatrix = worldMatrix * localMatrix.inverse();

  RenderMatrices renderMats;
  renderMats.localMatrix = localMatrix;
  renderMats.worldMatrix = worldMatrix;
  renderMats.cameraMatrix = cameraMatrix;
  renderMats.projMatrix = projMatrix;
  renderMats.screenMatrix = screenMatrix;
  return renderMats;
}

M4xX getScreenPoints(M4xX objSpacePts, RenderMatrices renderMats) {
  M4x4 localMatrixEigen = common::util::mMatrixToMatrix(renderMats.localMatrix);
  M4x4 worldMatrixEigen = common::util::mMatrixToMatrix(renderMats.worldMatrix);
  M4x4 cameraMatrixEigen =
      common::util::mMatrixToMatrix(renderMats.cameraMatrix);
  M4x4 projMatrixEigen = common::util::mMatrixToMatrix(renderMats.projMatrix);
  M4x4 screenMatrixEigen =
      common::util::mMatrixToMatrix(renderMats.screenMatrix);

  // Matrices for the pipeline
  M4xX objectToLocal = localMatrixEigen * objSpacePts;
  M4xX localToWorld = worldMatrixEigen * objectToLocal;
  M4xX worldToCam = cameraMatrixEigen * localToWorld;
  M4xX camToProj = projMatrixEigen * worldToCam;
  M4xX projToScreen = screenMatrixEigen * camToProj; // Might need to homogenize

  // Perspective Project in OpenGL (is this handled for orthographic?)
  // depths from zero to object's screen depth
  double alpha = projMatrixEigen(2, 2);
  double beta = projMatrixEigen(2, 3);

  VX camDepths = worldToCam.row(2);
  VX screenDepths = alpha * camDepths.array() + beta;
  camDepths *= -1.0; // This was negative before.

  // Get real obj screen points and stroke screen points depths?
  M4xX objScreenPts =
      projToScreen.array().rowwise() / camDepths.transpose().array();

  return objScreenPts;
}

} // namespace deformer
} // namespace abstract_squidget