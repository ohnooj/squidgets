#ifndef FACTORY_HPP
#define FACTORY_HPP

#include <maya/MDGModifier.h>
#include <maya/MDagPath.h>
#include <maya/MDoubleArray.h>
#include <maya/MObject.h>
#include <maya/MPointArray.h>

#include <maya/MFnNurbsCurve.h>

#include "common/util/maya.hpp"

namespace common {
namespace maya {
using namespace util;

MObject CreateNurbsCurve(MPointArray stroke, int deg,
                         MObject &parent = MObject::kNullObj);
} // namespace maya
} // namespace common

#endif
