#ifndef RENDER_MATRICES_HPP
#define RENDER_MATRICES_HPP

#include "abstract_squidget/types.hpp"

namespace abstract_squidget {
namespace deformer {

RenderMatrices getRenderMatrices(MObject geomObj);
M4xX getScreenPoints(M4xX objSpacePts, RenderMatrices renderMats);

} // namespace deformer
} // namespace abstract_squidget

#endif // RENDER_MATRICES_HPP