#ifndef ABSTRACT_DEFORMER_HPP
#define ABSTRACT_DEFORMER_HPP

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "abstract_squidget/types.hpp"
#include "common/logging/logger.hpp"

namespace abstract_squidget {
namespace deformer {
QueryClosestObjectResult queryPreselectObject(MString object,
                                              MPointArray strokePoints2D,
                                              MPointArray shadowCurve,
                                              EditFlags flags);

QueryClosestObjectResult
queryClosestObject(MStringArray candidateObjs, MPointArray curvePoints2D,
                   EditFlags flags,
                   std::shared_ptr<common::logging::Logger> logger);

DeformObjectResult transformObject(MString object, MPointArray curvePoints2D,
                                   RenderMatrices renderMats, EditFlags flags);

// Deform
DeformObjectResult deformObject(MString object, MPointArray curvePoints2D,
                                MPointArray shadowCurve,
                                RenderMatrices renderMats, EditFlags flags);

/**
 * @brief Performs a falloff on fnMesh.
 *
 * @param fnMesh FnMesh object to deform
 * @param handle_positions stroke vertices
 * @param start_stroke_i stroke index to start from (edge)
 * @param start_curve_i curve index to start from (edge)
 * @param interval_direction
 * @param renderMats
 */
void deformFalloff(MFnMesh &fnMesh, MPointArray handle_positions,
                   int start_stroke_i, int start_curve_i,
                   int interval_direction, RenderMatrices renderMats);

double sigmoid(double x);
double quintic(double x);

MMatrix screenTransformationMatrix(MVector t, double angle);

// Other
MVector livePoint(MString liveObject, MPoint mousePt, MVector offsetVec,
                  bool isCtrl);
MPointArray getObjSpaceCVs(MObject nurbsCurveObj, RenderMatrices renderMats);
} // namespace deformer
} // namespace abstract_squidget
#endif
