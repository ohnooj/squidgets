#ifndef TOON_RENDERER_HPP
#define TOON_RENDERER_HPP

#include <maya/M3dView.h>
#include <maya/MDagPath.h>
#include <maya/MDataHandle.h>
#include <maya/MGlobal.h>
#include <maya/MObject.h>
#include <maya/MObjectArray.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MPointArray.h>
#include <maya/MRenderLine.h>
#include <maya/MRenderLineArray.h>
#include <maya/MSelectionList.h>
#include <maya/MTransformationMatrix.h>
#include <maya/MVector.h>
#include <maya/MVectorArray.h>

#include <maya/MFnAttribute.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnPfxGeometry.h>
#include <maya/MFnTransform.h>

#include <common/maya/factory.hpp>

namespace common {
namespace maya {

class ToonRenderer {
public:
  ToonRenderer() {}
  ~ToonRenderer() {}

  MObjectArray renderCurves(MStringArray outputGeom, bool worldSpace);

private:
  void setupToonNode(MStringArray geomToRender);
  MObjectArray createToonCurves(bool worldSpace);

  MDagPath toonNodePath;
  MDagPath targetPath;
};

} // namespace maya
} // namespace common

#endif
