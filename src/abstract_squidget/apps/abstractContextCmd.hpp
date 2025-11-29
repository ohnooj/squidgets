#ifndef ABSTRACT_SQUIDGET_CONTEXT_CMD_HPP
#define ABSTRACT_SQUIDGET_CONTEXT_CMD_HPP

#include <maya/MPxContextCommand.h>

#include "abstract_squidget/apps/abstractContext.hpp"

namespace abstract_squidget {
namespace apps {
#define kStrokeWidthFlag "-sw"
#define kStrokeWidthFlagLong "-strokeWidth"
#define kStrokePointDistFlag "-d"
#define kStrokePointDistFlagLong "-distance"
#define kInputModeFlag "-im"
#define kInputModeFlagLong "-inputMode"
#define kTaskModeFlag "-tm"
#define kTaskModeFlagLong "-taskMode"
#define kCandidateObjectsFlag "-o"
#define kCandidateObjectsFlagLong "-objects"
#define kLoggerPathFlag "-lp"
#define kLoggerPathFlagLong "-loggerPath"
#define kSquidgetTypeFlag "-st"
#define kSquidgetTypeFlagLong "-squidgetType"

/**
 * Sets up CurveContext mouse tool.  Called from pluginMain.
 *
 * Can edit:
 * - AbstractContext::strokeWidth - Stroke width
 * - AbstractContext::strokePointDist - Stroke point distance
 * - AbstractContext::mode - study mode
 */
class AbstractContextCmd : public MPxContextCommand {
public:
  AbstractContextCmd();
  ~AbstractContextCmd();
  MPxContext *makeObj() override;
  static void *creator();

  MStatus doEditFlags() override;
  MStatus doQueryFlags() override;
  MStatus appendSyntax() override;

protected:
  AbstractContext *fAbstractContext;
};

} // namespace apps
} // namespace abstract_squidget

#endif