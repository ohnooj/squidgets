#include "bookmark_squidget/canvas/interface/element.hpp"

#include "common/util/math.hpp"

namespace bookmark_squidget {
namespace canvas {
namespace interface {

StepElement::StepElement() {
  renderConfig.color = STEP_ELEMENT_COLOR;
  renderConfig.strokeWidth = STEP_ELEMENT_STROKE_WIDTH;
  id = common::util::generateRandomID(4);
}

StepElement::~StepElement() {
  MGlobal::deleteNode(curve);
}

bool StepElement::isValid() {
  bool isNull = curve.isNull();
  return isNull;
}

} // namespace interface
} // namespace canvas
} // namespace bookmark_squidget