#ifndef ERROR_CALCULATION_HPP
#define ERROR_CALCULATION_HPP

#include <maya/MMatrix.h>
#include <maya/MPointArray.h>

#include "abstract_squidget/types.hpp"

namespace abstract_squidget {
namespace deformer {

// Error Calculation
double calculateShadowCurveError(RenderMatrices renderMats, MPointArray cvs,
                                 MPointArray curvePoints2D);
double matrixDeviation(MMatrix a, MMatrix b);
bool crossedOver(MPointArray curvePoints2D, MMatrix a, MMatrix b);
bool overRotate(MMatrix a, MMatrix b);

} // namespace deformer
} // namespace abstract_squidget

#endif // ERROR_CALCULATION_HPP