#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <string>

#include "maya/MColor.h"
#include "maya/MDagPath.h"
#include "maya/MGlobal.h"
#include "maya/MPoint.h"
#include "maya/MPointArray.h"

#include "maya/MFnDagNode.h"
#include "maya/MFnNurbsCurve.h"

#include "common/util/definitions.hpp"

namespace bookmark_squidget {
namespace canvas {
namespace interface {

namespace {
#define STEP_ELEMENT_COLOR MColor(1.0, 0.0, 0.0);
#define STEP_ELEMENT_STROKE_WIDTH 3.0;
} // namespace

struct RenderConfig {
  MColor color;
  double strokeWidth;
};

class Element {
  using Shape = MPointArray;

public:
  Shape shape;
  RenderConfig renderConfig;
};

class StepElement : public Element {
public:
  StepElement();
  ~StepElement();

  bool isValid();
  MObject curve;
  std::string id;
};

} // namespace interface
} // namespace canvas
} // namespace bookmark_squidget
#endif