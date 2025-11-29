#ifndef DEFORM_HPP
#define DEFORM_HPP

#include <Eigen/Dense>
#include <igl/AABB.h>

#include <maya/MDagPath.h>
#include <maya/MFnMesh.h>
#include <maya/MGlobal.h>
#include <maya/MPointArray.h>

#include "abstract_squidget/types.hpp"
#include "common/util/maya.hpp"

namespace abstract_squidget {
namespace deformer {
namespace method {
void findInterval(MPointArray meshVertices, MPointArray shadowCurve,
                  int &start_interval, int &end_interval, int &increment);
void convertPointsToScreenSpace(const MPointArray &objSpaceV,
                                const RenderMatrices &renderMats,
                                MPointArray &screenSpaceV);
void convertPointToScreenSpace(const MPoint &objSpaceV,
                               const RenderMatrices &renderMats,
                               MPoint &screenSpaceV);
Eigen::VectorXd parameterize_curve(const Eigen::MatrixX2d &curve);
Eigen::MatrixX2d convertMPointArrayToMatrixX2d(const MPointArray &mPointArray);
Eigen::VectorXd findPointOnCurve(const Eigen::MatrixXd &Y,
                                 const Eigen::VectorXd &Y_t, double t);

} // namespace method
} // namespace deformer
} // namespace abstract_squidget
#endif // DEFORM_HPP
