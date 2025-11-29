#ifndef BOOKMARKCONTEXTCMD_HPP
#define BOOKMARKCONTEXTCMD_HPP

#include <maya/MArgParser.h>
#include <maya/MPxContextCommand.h>
#include <maya/MSyntax.h>

#include "bookmark_squidget/apps/bookmarkContext.hpp"

namespace bookmark_squidget {
namespace apps {

#define kStrokeWidthFlag "-sw"
#define kStrokeWidthFlagLong "-strokeWidth"
#define kStrokePointDistFlag "-d"
#define kStrokePointDistFlagLong "-distance"
// #define kInputModeFlag "-im"
// #define kInputModeFlagLong "-inputMode"
// #define kTaskModeFlag "-tm"
// #define kTaskModeFlagLong "-taskMode"
// #define kCandidateObjectsFlag "-o"
// #define kCandidateObjectsFlagLong "-objects"
// #define kLoggerPathFlag "-lp"
// #define kLoggerPathFlagLong "-loggerPath"

/**
 * Sets up BookmarkContext mouse tool.  Called from pluginMain.
 *
 * Can edit:
 * - BookmarkContext::strokeWidth
 * - BookmarkContext::strokePointDist
 *
 */
class BookmarkContextCmd : public MPxContextCommand {
public:
  BookmarkContextCmd();
  ~BookmarkContextCmd();
  MPxContext *makeObj() override;
  static void *creator();

  MStatus doEditFlags() override;
  MStatus doQueryFlags() override;
  MStatus appendSyntax() override;

protected:
  BookmarkContext *fBookmarkContext;
};

} // namespace apps
} // namespace bookmark_squidget

#endif
