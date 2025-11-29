#include "bookmark_squidget/apps/bookmarkSquidgetTool.hpp"

namespace bookmark_squidget {
namespace apps {
BookmarkSquidgetTool::CurveDistanceItem BookmarkSquidgetTool::toDeformItem =
    BookmarkSquidgetTool::CurveDistanceItem();

BookmarkSquidgetTool::BookmarkSquidgetTool() {
  setCommandString("BookmarkSquidgetToolCmd");
}

BookmarkSquidgetTool::~BookmarkSquidgetTool() {}

void *BookmarkSquidgetTool::creator() { return new BookmarkSquidgetTool; }

MStatus BookmarkSquidgetTool::doIt(const MArgList &args) {
  MStatus status = parseArgs(args);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  predoIt();
  return redoIt();
}

MStatus BookmarkSquidgetTool::predoIt() {
  MStatus status = MS::kSuccess;
  CHECK_MSTATUS_AND_RETURN_IT(status);

  cout << "BookmarkSquidgetTool::predoIt() " << toDeformItem.score << endl;
  if (toDeformItem.score == -1) {
    queryClosestObjectOnly();
  }

  return status;
}

MStatus BookmarkSquidgetTool::cancel() { return MStatus::kSuccess; }

MStatus BookmarkSquidgetTool::parseArgs(const MArgList &args) {
  MStatus status;
  MArgDatabase argData(syntax(), args);

  return status;
}

MStatus BookmarkSquidgetTool::redoIt() {
  MStatus status;
  // generate deformation map
  deformClosestObject();

  cout << "after deform" << endl;
  toDeformItem = CurveDistanceItem(); // Reset
  // common::util::log::implicitSet(controlPath.fullPathName().asChar());
  return status;
}

MStatus BookmarkSquidgetTool::undoIt() {
  MStatus status;
  // common::util::log::log("impl-<UNDO>");

  MFnTransform fnTransform(controlPath, &status);
  fnTransform.set(oldMatrix);

  controlPath = MDagPath();
  return status;
}

MStatus BookmarkSquidgetTool::finalize() {
  MArgList command;
  command.addArg(commandString());
  return MPxToolCommand::doFinalize(command);
}

MSyntax BookmarkSquidgetTool::newSyntax() {
  MSyntax syntax;
  return syntax;
}

// =============================================================================
//                            Deform Object
// =============================================================================
// void BookmarkSquidgetTool::queryClosestObjectOnly(MPointArray strokePoints2D)
// {
void BookmarkSquidgetTool::queryClosestObjectOnly() {
  /**
   * @brief Query closest object to the 2D stroke points.
   * Closest object is defined by the following:
   * - Distance from stroke points to object's shadow curve.
   * - Shape difference of object's shadow curve to stroke points.
   */
  MStringArray outputGeoms = getGeomsAroundStroke(); // Select closest Geometry

  std::priority_queue<CurveDistanceItem, std::vector<CurveDistanceItem>,
                      CompareCurveDistanceItem>
      pq;

  for (MString outputGeom : outputGeoms) { // For each output geometry
    MObject geomObj = common::util::getObjectByName(outputGeom);
    MFnDagNode geomFn(geomObj);
    MDagPath geomPath;
    geomFn.getPath(geomPath);

    MString inputGeoms[] = {outputGeom};
    MObjectArray curves = getImplicitSquidgetCurves(
        MStringArray(inputGeoms, 1), false); // Generate toon curves
    // curves = common::util::combineAndFilterCurves(curves);

    for (MObject curveObj : curves) { // Eigen Matrices
      RenderMatrices renderMats = getRenderMatrices(geomObj);
      MMatrix initiallocalMatrix = renderMats.localMatrix;

      MFnDagNode fnDagNode(curveObj);
      MDagPath curvePath;
      fnDagNode.getPath(curvePath);
      MFnNurbsCurve curveFn(curvePath);

      MPointArray curveWorldCVs;
      curveFn.getCVs(curveWorldCVs);

      MPointArray curveObjCVs;
      for (MPoint p : curveWorldCVs) { // Convert CVs to object space to match
                                       // object matrix
        MPoint objp =
            (renderMats.worldMatrix * renderMats.localMatrix).inverse() * p;
        curveObjCVs.append(objp);
      }

      MPointArray objSpaceShadowCurve =
          findShadowCurve( // Find Object Space Shadow Curve
              penControlPoints, curveObjCVs, renderMats, 0);

      // A S P C W L O = S P C W B L O
      // A X L O = X B L O
      // A X = X B
      // X^-1 A X = B
      // M4x4 A = leastSquares(curveObjCVs, penControlPoints, renderMats); //
      // Least squares at screen space
      M4x4 A = leastSquares(objSpaceShadowCurve, penControlPoints,
                            renderMats); // Least squares at screen space
      M4x4 X = common::util::mMatrixToMatrix(renderMats.screenMatrix) *
               common::util::mMatrixToMatrix(renderMats.projMatrix) *
               common::util::mMatrixToMatrix(renderMats.cameraMatrix) *
               common::util::mMatrixToMatrix(renderMats.worldMatrix);
      M4x4 B = X.inverse() * A * X;

      renderMats.screenTransformMatrix = common::util::matrixToMMatrix(
          A); // How to update the screen transform matrix
      MMatrix newLocalMatrix =
          common::util::matrixToMMatrix(B) *
          renderMats.localMatrix; // Update local matrix with new screen
                                  // transform matrix

      // Project into screen space
      // I don't use worldToScreen because i don't remember why, but this is
      // more accurate.
      double squared_norm = calculateShadowCurveError(renderMats, curveObjCVs);
      double matrix_distance =
          matrixDeviation(initiallocalMatrix, newLocalMatrix);

      // pq.emplace(squared_norm, matrix_distance, renderMats, geomPath,
      // objSpaceShadowCurve);
      pq.emplace(squared_norm, matrix_distance, renderMats, geomPath,
                 curveObjCVs);
    }

    for (MObject curveObj : curves) { // Delete intermediary toon curves.
      common::util::deleteCurve(curveObj);
    }
  }

  if (pq.size() == 0) { // cout << "pq is empty" << endl;
    toDeformItem = CurveDistanceItem();
    return;
  } else {
    toDeformItem = pq.top();
    MGlobal::selectByName(toDeformItem.geomPath.fullPathName(),
                          MGlobal::kReplaceList);
  }
}

void BookmarkSquidgetTool::deformClosestObject() {
  cout << "deformClosestObject" << endl;
  RenderMatrices renderMats = getRenderMatrices(toDeformItem.geomPath.node());
  MPointArray curveObjCVs = toDeformItem.shadowCurve; // Not shadow curve rly

  double noise = 0; // Iterative Least Squares with regularization
  int runs = 5;
  for (int i = 0; i < runs; ++i) {
    double runNoise = noise * (1 - (double)i / runs);
    MPointArray objSpaceShadowCurve =
        findShadowCurve(penControlPoints, curveObjCVs, renderMats,
                        runNoise); // Find Object Space Shadow Curve

    // common::maya::CreateNurbsCurve(objSpaceShadowCurve, 3);

    // A S P C W L O = S P C W B L O
    // A X L O = X B L O
    // A X = X B
    // X^-1 A X = B
    M4x4 A = leastSquares(objSpaceShadowCurve, penControlPoints,
                          renderMats); // Least squares at screen space
    M4x4 X = common::util::mMatrixToMatrix(renderMats.screenMatrix) *
             common::util::mMatrixToMatrix(renderMats.projMatrix) *
             common::util::mMatrixToMatrix(renderMats.cameraMatrix) *
             common::util::mMatrixToMatrix(renderMats.worldMatrix);
    M4x4 B = X.inverse() * A * X;

    renderMats.screenTransformMatrix = common::util::matrixToMMatrix(A);
    renderMats.localMatrix =
        common::util::matrixToMMatrix(B) * renderMats.localMatrix;
  }
  // common::util::deleteCurve(temp);
  cout << "loops" << endl;

  MFnTransform fnMinTransform;
  fnMinTransform.setObject(toDeformItem.geomPath);
  fnMinTransform.getPath(controlPath);
  controlMatrices = renderMats;

  oldMatrix =
      fnMinTransform.transformationMatrix(); // Save old transform for undo
  fnMinTransform.set(
      MTransformationMatrix(controlMatrices.localMatrix.transpose())
          .asMatrix());
  // fnMinTransform.setScale(scale); // Maintain scale?
  double sheer[3] = {0, 0, 0}; // Maintain Sheer
  fnMinTransform.setShear(sheer);

  MStatus s;
  controlPoint = fnMinTransform.getTranslation(MSpace::kPostTransform,
                                               &s); // Perform transformation.
  MPoint controlScreenPt =
      (controlMatrices.screenMatrix * controlMatrices.projMatrix *
       controlMatrices.cameraMatrix * controlMatrices.worldMatrix) *
      controlPoint;
  controlScreenPt.cartesianize();

  oldDragDifference =
      penControlPoints[penControlPoints.length() - 1] -
      controlScreenPt; // Save old drag difference for screen space holding.
  oldDragDifference.z = 0;
}

MStringArray BookmarkSquidgetTool::getGeomsAroundStroke() {
  MStringArray outputGeoms;
  if (mode == 0) {
    outputGeoms.append("pTorus1");
    outputGeoms.append("pSphere1");
    outputGeoms.append("pCube1");
    outputGeoms.append("pCylinder1");
    outputGeoms.append("pSuperShape1");
  } else if (mode == 1) {
    outputGeoms.append("moon");
    outputGeoms.append("trash");
    outputGeoms.append("trash1");
    outputGeoms.append("trash2");
    outputGeoms.append("lamp1");
    outputGeoms.append("lamp2");
    outputGeoms.append("lamp3");
    outputGeoms.append("Bench1");
    outputGeoms.append("Bench2");
    outputGeoms.append("Bench3");
    outputGeoms.append("Bench4");
    outputGeoms.append("Bench5");
  } else if (mode == 2) {
    // outputGeoms.append("lou:FKShoulder_R");
    // outputGeoms.append("lou:FKElbow_R");
    // outputGeoms.append("lou:FKWrist_R");
    // outputGeoms.append("lou:FKShoulder_L");
    // outputGeoms.append("lou:FKElbow_L");
    // outputGeoms.append("lou:FKWrist_L");
    // outputGeoms.append("lou:FKRoot_M");
    // outputGeoms.append("lou:FKSpine1_M");
    // outputGeoms.append("lou:FKChest_M");
    // outputGeoms.append("lou:FKNeck_M");
    // outputGeoms.append("lou:FKHead_M");
    // outputGeoms.append("lou:Main");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_root_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_spine_ik_hip_ctrl");
    // outputGeoms.append("squirrels_shot_ready_anim:sqr_spine_curvature_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_spine_ik_shoulder_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_lf_ear_2_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_rt_ear_2_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_tail_1_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_tail_2_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_tail_3_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_tail_4_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_tail_5_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_tail_6_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_tail_7_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_lower_jaw_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_upper_jaw_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_rt_mouth_corner_ctrl");
    outputGeoms.append("squirrels_shot_ready_anim:sqr_lf_mouth_corner_ctrl");
  }
  return outputGeoms;
}

MObjectArray
BookmarkSquidgetTool::getImplicitSquidgetCurves(MStringArray inputGeom,
                                                bool worldSpace) {
  MObjectArray curves;

  MSelectionList oldList;
  MGlobal::getActiveSelectionList(oldList);

  if (mode == 0 || mode == 1) { // For regular shapes, use toon renderer
    common::maya::ToonRenderer r;
    curves = r.renderCurves(inputGeom, worldSpace);
  }

  if (mode ==
      2) { // For pre-authored implicit squidgets, use curves from scene.
    for (MString g : inputGeom) {
      MGlobal::selectByName(g, MGlobal::kReplaceList);

      MSelectionList sel_list;
      MGlobal::getActiveSelectionList(sel_list);

      MDagPath geomPath;
      MFnNurbsCurve curveFn;
      MPointArray cvs;
      sel_list.getDagPath(0, geomPath);
      geomPath.extendToShapeDirectlyBelow(0);
      curveFn.setObject(geomPath);
      curveFn.getCVs(cvs, MSpace::kWorld);

      MObject copyCurve =
          common::maya::CreateNurbsCurve(cvs, 3); // Create fake implicit curves
      curves.append(copyCurve);

      //  MFnNurbsCurve ocopyCurveFn;
      // MObject copyCurve = copyCurveFn.copy(geomPath.node());
      // copyCurveFn.setCVs(cvs);
      // MDagPath copyPath;
      // copyCurveFn.getPath(copyPath);
      // curves.append(copyPath.node());
    }
  }

  // MGlobal::setActiveSelectionList(oldList);
  return curves;
}

BookmarkSquidgetTool::RenderMatrices
BookmarkSquidgetTool::getRenderMatrices(MObject geomObj) {
  // cout << "BookmarkSquidgetTool::getRenderMatrices()" << endl;
  MFnDagNode fnDagNode(geomObj);
  MDagPath geomPath;
  fnDagNode.getPath(geomPath);
  MFnTransform fnTransform(geomPath);
  // cout << fnTransform.fullPathName() << endl;

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
  // worldMat.setRotationQuaternion(worldRot.x, worldRot.y, worldRot.z,
  // worldRot.w, MSpace::kWorld);
  worldMat.setTranslation(worldTrans, MSpace::kWorld);

  // World = Parent * Local
  // World * Local^-1 = Parent
  MMatrix worldMatrix =
      worldMat.asMatrix().transpose(); // This is actually parent matrix
  worldMatrix = worldMatrix * localMatrix.inverse();

  // cout << "worldMatrix localMatrix" << endl;
  // cout << worldMatrix << endl;
  // cout << localMatrix << endl;

  RenderMatrices renderMats;
  renderMats.localMatrix = localMatrix;
  renderMats.worldMatrix = worldMatrix;
  renderMats.cameraMatrix = cameraMatrix;
  renderMats.projMatrix = projMatrix;
  renderMats.screenMatrix = screenMatrix;
  return renderMats;
}

MPointArray BookmarkSquidgetTool::findShadowCurve(MPointArray strokePoints2D,
                                                  MPointArray cvs,
                                                  RenderMatrices renderMats,
                                                  double noise) {
  /**
   * @brief pojects strokePoints2D onto cvs using curveObj MFnCurve parameters
   * in screen space.  Returns the projected points in object space.
   *
   * DOESN"T RETURN PROJECTION.  RETURNS SAMPLED OUTLINE NOW.
   */

  MObject curveObj = common::maya::CreateNurbsCurve(cvs, 0);
  MFnDagNode fnDagNode(curveObj);
  MDagPath curvePath;
  fnDagNode.getPath(curvePath);
  MFnNurbsCurve curveFn(curvePath);

  // Sample curve to have more control points on the curve.
  // cout << "cvs:" << endl << cvs << endl;
  MPointArray sampled_cvs;
  int samples = 100;
  for (int i = 0; i < samples; ++i) {
    // double param = (double)i / (samples - 1);
    double param = cvs.length() * (double)i / (samples - 1);
    MPoint p;
    curveFn.getPointAtParam(param, p);
    sampled_cvs.append(p);
  }
  // cout << "sampled_cvs:" << endl << sampled_cvs << endl;

  MPointArray worldCVs;
  MPointArray screenCVs;
  for (MPoint p : sampled_cvs) {
    MPoint worldp = renderMats.worldMatrix * renderMats.localMatrix * p;
    worldCVs.append(worldp);

    MPoint screenp = renderMats.screenMatrix * renderMats.projMatrix *
                     renderMats.cameraMatrix * worldp;
    screenp.cartesianize();
    screenCVs.append(screenp);
  }

  // Project screen points onto screenCVs to map stroke points to curve.
  // But rather we get the parameter of the closest point on the curve to the
  // stroke point. Such that we can find the corresponding worldCVs point using
  // the parameter. Add random noise to point to account for symmetry.  Escape
  // from equillibrium?
  curveFn.setCVs(screenCVs);
  MDoubleArray shadowParam;
  for (MPoint p : strokePoints2D) {
    double param;
    MPoint noisyP =
        p + noise * (MVector(rand() % 100, rand() % 100, rand() % 100) / 100.);
    MPoint shadowPoint = curveFn.closestPoint(noisyP, &param);
    shadowParam.append(param);
  }

  // Find shadow curve in world space and then convert to object space.
  curveFn.setCVs(worldCVs); // Reset original cvs
  MPointArray objSpaceShadowCurve;
  for (double param : shadowParam) {
    MPoint shadowCVPoint;
    curveFn.getPointAtParam(param, shadowCVPoint);
    MPoint objp = (renderMats.worldMatrix * renderMats.localMatrix).inverse() *
                  shadowCVPoint;
    objSpaceShadowCurve.append(objp);
  }

  common::util::deleteCurve(curveObj);
  // return objSpaceShadowCurve;
  return sampled_cvs;
}

// =============================================================================
//                            Least Squares
// =============================================================================
Eigen::Matrix3d BookmarkSquidgetTool::leastSquaresTemp(M3xX x0, M3xX b0) {
  // Calculate Polar Decomposition for transformation matrix
  M2xX x = x0.block(0, 0, 2, x0.cols());
  M2xX b = b0.block(0, 0, 2, b0.cols());

  V2 xAvgPoint = V2::Zero();
  V2 bAvgPoint = V2::Zero();

  for (int i = 0; i < x.cols(); ++i) { // cout << "x and b avg" << endl;
    xAvgPoint += x.col(i);
    bAvgPoint += b.col(i);
  }
  xAvgPoint /= x.cols();
  bAvgPoint /= b.cols();

  M2xX xLocal = x.colwise() - xAvgPoint; // cout  << "x and b local" << endl;
  M2xX bLocal = b.colwise() - bAvgPoint;

  M2x2 Apq = bLocal * xLocal.transpose(); // cout << "Apq and Aqq" << endl;
  M2x2 Aqq = xLocal * xLocal.transpose();

  Eigen::JacobiSVD<M2x2> svd(Apq, Eigen::ComputeFullU | Eigen::ComputeFullV);
  M2x2 U = svd.matrixU();
  M2x2 V = svd.matrixV();
  M2x2 D = svd.singularValues().asDiagonal();

  M2x2 polar0 = U * V.transpose();
  V2 trans = bAvgPoint - polar0 * xAvgPoint;

  M3x3 ret = M3x3::Zero();
  ret.block(0, 0, 2, 2) = polar0; // Set top left 2x2 of ret to polar0
  ret.block(0, 2, 2, 1) = trans;  // Set top right of ret to translation
  ret(2, 2) = 1;                  // Set bottom right 1x1 of ret to 1

  return ret;
}

Eigen::Matrix4d BookmarkSquidgetTool::leastSquares(MPointArray initialPoints,
                                                   MPointArray targetPoints,
                                                   RenderMatrices renderMats) {
  // cout << "BookmarkSquidgetTool::leastSquares()" << endl;
  M4x4 localMatrixEigen = common::util::mMatrixToMatrix(renderMats.localMatrix);
  M4x4 worldMatrixEigen = common::util::mMatrixToMatrix(renderMats.worldMatrix);
  M4x4 cameraMatrixEigen =
      common::util::mMatrixToMatrix(renderMats.cameraMatrix);
  M4x4 projMatrixEigen = common::util::mMatrixToMatrix(renderMats.projMatrix);
  M4x4 screenMatrixEigen =
      common::util::mMatrixToMatrix(renderMats.screenMatrix);

  // cout << "o and v" << endl;
  M4xX o = common::util::pointArrayToMatrix(initialPoints); // Solve Ao = v
  M4xX v = common::util::pointArrayToMatrix(targetPoints);

  M4xX objectToLocal = localMatrixEigen * o; // Object to screen
  M4xX localToWorld = worldMatrixEigen * objectToLocal;
  M4xX worldToCam = cameraMatrixEigen * localToWorld;
  M4xX camToProj = projMatrixEigen * worldToCam;
  M4xX projToScreen = screenMatrixEigen * camToProj; // Might need to homogenize

  // cout << "alpha and beta" << endl;
  double alpha = projMatrixEigen(
      2, 2); // Perspective Project in OpenGL       // Set PenScreenPoints
             // depths from zero to object's screen depth
  double beta = projMatrixEigen(2, 3);

  VX camDepths = worldToCam.row(2);
  VX screenDepths = alpha * camDepths.array() + beta;
  camDepths *= -1.0; // This was negative before.

  M4xX objScreenPts =
      projToScreen.array().rowwise() /
      camDepths.transpose().array(); // Get real obj screen points and stroke
                                     // screen points depths
  // cout << "objScreenPts and v" << endl;
  // v.row(2) = screenDepths.array() / camDepths.array(); // Set depth for
  // screen points since this is unknown.

  M3xX objScreen2d = objScreenPts.block(
      0, 0, 3, objScreenPts.cols()); // Convert 3d points to 2d points.
  M3xX v2d = v.block(0, 0, 3, v.cols());

  objScreen2d.row(2) = VX::Ones(objScreen2d.cols());
  v2d.row(2) = VX::Ones(v2d.cols());

  // cout << "objScreen2d and v2d" << endl;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> distribution(
      0.0, 0.1); // Mean=0, Standard Deviation=0.01
  // Add small noise to each element in the matrix
  for (int i = 0; i < objScreen2d.rows() - 1; ++i) {
    for (int j = 0; j < objScreen2d.cols(); ++j) {
      objScreen2d(i, j) += distribution(gen);
    }
  }

  for (int i = 0; i < v2d.rows() - 1; ++i) {
    for (int j = 0; j < v2d.cols(); ++j) {
      v2d(i, j) += distribution(gen);
    }
  }

  // cout << "coordinating" << endl;
  Eigen::MatrixX2d coordinates_1 =
      objScreen2d.block(0, 0, 2, objScreen2d.cols()).transpose();
  Eigen::Vector2d centroid_1 = Eigen::Vector2d::Zero();
  for (int i = 0; i < coordinates_1.rows(); ++i) {
    centroid_1 += coordinates_1.row(i);
  }
  centroid_1 /= coordinates_1.rows();
  // cout << "coordinates_1: " << endl;
  // for (int i = 0; i < coordinates_1.rows(); ++i) {
  //   cout << "[" << coordinates_1(i, 0) << ", " << coordinates_1(i, 1) << "],"
  //   << endl;
  // }

  // Initialize the second matrix with 34 rows and 2 columns
  Eigen::MatrixX2d coordinates_2 = v2d.block(0, 0, 2, v2d.cols()).transpose();
  Eigen::Vector2d centroid_2 = Eigen::Vector2d::Zero();
  for (int i = 0; i < coordinates_2.rows(); ++i) {
    centroid_2 += coordinates_2.row(i);
  }
  centroid_2 /= coordinates_2.rows();
  // cout << "coordinates_2: " << endl;
  // for (int i = 0; i < coordinates_2.rows(); ++i) {
  //   cout << "[" << coordinates_2(i, 0) << ", " << coordinates_2(i, 1) << "],"
  //   << endl;
  // }

  double dist = common::algo::euclidean_distance(centroid_1, centroid_2);
  // cout << "dist: " << dist << endl;
  // cout << "objScreen2d: " << objScreen2d.size() << endl;
  // cout << "v2d: " << v2d.size() << endl;
  std::pair<std::vector<Eigen::Matrix3d>, Eigen::MatrixX2d> icpResult =
      common::algo::icp(coordinates_1, coordinates_2, 20, 2 * dist, 1e-3, 1e-4,
                        10, false);

  M4x4 icpAll = M4x4::Identity(); // Homogeneous Matrix Convert from 3x3 to 4x4
  for (size_t i = 0; i < icpResult.first.size(); ++i) {
    M4x4 icp = M4x4::Identity();
    icp.block(0, 0, 2, 2) = icpResult.first[i].block(0, 0, 2, 2);
    icp.block(0, 3, 2, 1) = icpResult.first[i].block(0, 2, 2, 1);
    // cout << i << ": " << endl << icp << endl;
    icpAll = icp * icpAll;
  }

  // cout << "icpAll: " << endl << icpAll << endl;
  M4x4 icpAllInv = icpAll.inverse();
  // cout << "icpAllInv: " << endl << icpAllInv << endl;

  // M3x3 screenA = leastSquaresTemp(objScreen2d, v2d); // Least Squares
  // calculation for screen space cout << "screenA: " << endl << screenA <<
  // endl; M4x4 screenA4 = M4x4::Identity(); // Homogeneous Matrix Convert from
  // 3x3 to 4x4 screenA4.block(0, 0, 2, 2) = screenA.block(0, 0, 2, 2);
  // screenA4.block(0, 3, 2, 1) = screenA.block(0, 2, 2, 1);

  // M4x4 screenA4 = M4x4::Identity(); // Homogeneous Matrix Convert from 3x3 to
  // 4x4
  M4x4 screenA4 = icpAllInv; // Homogeneous Matrix Convert from 3x3 to 4x4
  // cout << "screenA4: " << endl << screenA4 << endl;
  // A S P C W = S P C B W
  // A X W  = X B W
  // A X = X B
  // X^-1 A X = B
  // M4x4 A = screenA4;
  // M4x4 X = screenMatrixEigen * projMatrixEigen * cameraMatrixEigen;
  // M4x4 B = X.inverse() * A * X;
  // // //cout << "A: " << endl << A << endl;
  // // //cout << "B: " << endl << B << endl;
  // return B;

  return screenA4;
}

double
BookmarkSquidgetTool::calculateShadowCurveError(RenderMatrices renderMats,
                                                MPointArray cvs) {
  MPointArray screenCVs;
  for (MPoint p : cvs) {
    MPoint worldp = renderMats.screenMatrix * renderMats.projMatrix *
                    renderMats.cameraMatrix * renderMats.worldMatrix * p;
    worldp.cartesianize();
    screenCVs.append(worldp);
  }

  // Find closest point of screenCVs on curve to strokePoints
  // Find params of closest points.  Map params to curve.
  MObject curveObj = common::maya::CreateNurbsCurve(cvs, 0);
  MFnDagNode fnDagNode(curveObj);
  MDagPath curvePath;
  fnDagNode.getPath(curvePath);
  MFnNurbsCurve curveFn(curvePath);
  // curveFn.setCVs(screenCVs);
  MPointArray shadowCurve;
  for (MPoint p : penControlPoints) {
    MPoint shadowPoint = curveFn.closestPoint(p);
    shadowCurve.append(shadowPoint);
  }
  common::util::deleteCurve(curveObj);

  double squared_norm = 0;
  for (unsigned int i = 0; i < shadowCurve.length(); ++i) {
    squared_norm += pow((shadowCurve[i] - penControlPoints[i]).length(), 2);
  }
  squared_norm /= shadowCurve.length();
  return squared_norm;
}

double BookmarkSquidgetTool::matrixDeviation(MMatrix a, MMatrix b) {
  MVector aT(a[0][3], a[1][3], a[2][3]);
  MVector bT(b[0][3], b[1][3], b[2][3]);

  return pow((aT - bT).length(), 2);
}

bool BookmarkSquidgetTool::crossedOver(MMatrix a, MMatrix b) {
  bool crossesOver = false;
  MPoint aT(a[0][3], a[1][3], a[2][3]);
  MPoint bT(b[0][3], b[1][3], b[2][3]);

  M3dView view = M3dView::active3dView();
  short ax, ay, bx, by;
  view.worldToView(aT, ax, ay);
  view.worldToView(bT, bx, by);
  MPoint aScreen(ax, ay, 0);
  MPoint bScreen(bx, by, 0);

  MPoint strokeMidpoint(0, 0, 0);
  for (MPoint p : penControlPoints) {
    strokeMidpoint += p;
  }
  strokeMidpoint = strokeMidpoint / penControlPoints.length();

  MVector AP = aScreen - strokeMidpoint;
  MVector BP = bScreen - strokeMidpoint;

  double angle = AP.angle(BP);
  // cout << "angle: " << angle << " " << aScreen << bScreen << strokeMidpoint
  // << endl;
  return angle > M_PI / 2.;
}

bool BookmarkSquidgetTool::overRotate(MMatrix a, MMatrix b) {
  // checks if rotation of a to b is greater than 90 degrees
  MTransformationMatrix aTrans(a);
  MTransformationMatrix bTrans(b);

  MVector aRot = aTrans.eulerRotation().asVector();
  MVector bRot = bTrans.eulerRotation().asVector();

  double angle = aRot.angle(bRot);
  // //cout << "angle: " << angle << endl;
  return angle > M_PI / 2.;
}

// =============================================================================
//                            Other
// =============================================================================
void BookmarkSquidgetTool::ReadPoint(MPoint mousePt, MVector offsetVec,
                                     bool isCtrl) {
  // Assume controlObj exists
  if (controlPath.isValid()) {
    MFnTransform fnTransform(controlPath);

    MStatus s;
    // //cout << fnTransform.fullPathName() << endl;
    MPoint controlPoint = fnTransform.getTranslation(MSpace::kWorld, &s);
    MPoint controlScreenPt =
        (controlMatrices.screenMatrix * controlMatrices.projMatrix *
         controlMatrices.cameraMatrix) *
        controlPoint;

    M4x4 screenATrans = M4x4::Identity();
    M4x4 screenARot = M4x4::Identity();

    M4x4 T = M4x4::Identity();
    M4x4 R = M4x4::Identity();
    if (isCtrl) {
      // common::util::log::implicitHold(controlPath.fullPathName().asChar(),
      //                                 "Modifier");
      MVector difference = mousePt - controlScreenPt;
      // //cout << "Mouse: " << mousePt << " "
      // << controlScreenPt << " -- " <<
      // difference << endl;
      difference.z = 0;
      // Check if difference is rotated
      // clockwise to oldDragDifference
      double theta = difference.angle(oldDragDifference);
      bool isClockwise = difference.x * oldDragDifference.y -
                             difference.y * oldDragDifference.x >
                         0;
      if (isClockwise) {
        theta = -theta;
      }

      // //cout << "theta: " << theta << " " <<
      // difference << " | " <<
      // oldDragDifference << endl;
      oldDragDifference = difference;

      screenARot(0, 0) = cos(theta);
      screenARot(0, 1) = -sin(theta);
      screenARot(1, 0) = sin(theta);
      screenARot(1, 1) = cos(theta);
      // cout << screenARot << endl;

      // A S' S P C W = S' S P C (T W R)
      // A S' S P C W = S' S P C (T W R)
      // A X = X R
      // X^-1 A X = R
      M4x4 A = screenARot;
      M4x4 screenTranslate = M4x4::Identity();
      screenTranslate(0, 3) = -controlMatrices.screenMatrix(0, 3);
      screenTranslate(1, 3) = -controlMatrices.screenMatrix(1, 3);

      M4x4 localTranslate = M4x4::Identity();
      localTranslate(0, 3) = -controlMatrices.localMatrix(0, 3);
      localTranslate(1, 3) = -controlMatrices.localMatrix(1, 3);
      localTranslate(2, 3) = -controlMatrices.localMatrix(2, 3);

      M4x4 X = screenTranslate *
               common::util::mMatrixToMatrix(controlMatrices.screenMatrix) *
               common::util::mMatrixToMatrix(controlMatrices.projMatrix) *
               common::util::mMatrixToMatrix(controlMatrices.cameraMatrix) *
               common::util::mMatrixToMatrix(controlMatrices.worldMatrix) *
               localTranslate;
      R = X.inverse() * A * X;

      // A S' S P C = S' S P C -T R T
      // X =
      // common::util::mMatrixToMatrix(controlMatrices.screenMatrix)
      // *
      //       common::util::mMatrixToMatrix(controlMatrices.projMatrix)
      //       *
      //       common::util::mMatrixToMatrix(controlMatrices.cameraMatrix);
      // R = X.inverse() * A * X;
      MVector oldTrans;
      oldTrans.x = controlMatrices.localMatrix(0, 3);
      oldTrans.y = controlMatrices.localMatrix(1, 3);
      oldTrans.z = controlMatrices.localMatrix(2, 3);

      controlMatrices.localMatrix =
          common::util::matrixToMMatrix(R) * controlMatrices.localMatrix;
      controlMatrices.localMatrix(0, 3) = oldTrans.x;
      controlMatrices.localMatrix(1, 3) = oldTrans.y;
      controlMatrices.localMatrix(2, 3) = oldTrans.z;

    } else {
      // common::util::log::implicitHold(controlPath.fullPathName().asChar(),
      // ""); Translate only by mouseObjectOffsetPoint and currentObjectPoint
      // MPoint offsetMousePt = mousePt + offsetVec;
      MPoint offsetMousePt = mousePt - oldDragDifference;
      MVector difference = offsetMousePt - controlScreenPt;
      screenATrans(0, 3) = difference.x;
      screenATrans(1, 3) = difference.y;

      // A S P C W L = S P C W (T L R)
      // A S P C W L = S P C W T L
      // A X = X T
      // X^-1 A X = T
      M4x4 A = screenATrans;
      M4x4 X = common::util::mMatrixToMatrix(controlMatrices.screenMatrix) *
               common::util::mMatrixToMatrix(controlMatrices.projMatrix) *
               common::util::mMatrixToMatrix(controlMatrices.cameraMatrix) *
               common::util::mMatrixToMatrix(controlMatrices.worldMatrix);
      T = X.inverse() * A * X;
      controlMatrices.localMatrix =
          common::util::matrixToMMatrix(T) * controlMatrices.localMatrix;
    }

    fnTransform.set(
        MTransformationMatrix(controlMatrices.localMatrix.transpose())
            .asMatrix());
    // double scale[3] = {1, 1, 1};
    // fnTransform.setScale(scale);
    double sheer[3] = {0, 0, 0};
    fnTransform.setShear(sheer);
  }
}

MStatus BookmarkSquidgetTool::getControlPt(MPoint &controlPt) {
  controlPt = controlPoint;
  return MStatus::kSuccess;
}

void BookmarkSquidgetTool::setPenControlPoints(MPointArray arr) {
  penControlPoints = arr;
}

} // namespace apps
} // namespace bookmark_squidget
