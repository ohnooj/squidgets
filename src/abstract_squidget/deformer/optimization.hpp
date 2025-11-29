#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

#include <Eigen/Dense>
#include <maya/MPointArray.h>

#include "abstract_squidget/types.hpp"

namespace abstract_squidget {
namespace deformer {

RenderMatrices optimizeNewMatrices(MString object, MPointArray strokePoints2D,
                                   MPointArray shadowCurve, EditFlags flags);

RenderMatrices optimizeNewMatrices(MPointArray squidgetObjCVs,
                                   MPointArray curvePoints2D,
                                   RenderMatrices renderMats, EditFlags flags);

// Least Squares
Eigen::Matrix4d leastSquares(MPointArray initialPoints,
                             MPointArray targetPoints,
                             RenderMatrices renderMats, EditFlags);

Eigen::Matrix4d optimizeTranslation(M3xX v2d, M3xX objScreen2d,
                                    EditFlags flags);

} // namespace deformer
} // namespace abstract_squidget

#endif // OPTIMIZATION_HPP