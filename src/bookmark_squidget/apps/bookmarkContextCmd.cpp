#include "bookmark_squidget/apps/bookmarkContextCmd.hpp"

namespace bookmark_squidget {
namespace apps {

BookmarkContextCmd::BookmarkContextCmd() {
  // cout << "BookmarkContextCmd()" << endl;
}

BookmarkContextCmd::~BookmarkContextCmd() {
  // cout << "~BookmarkContextCmd()" << endl;
  delete fBookmarkContext;
}

MStatus BookmarkContextCmd::doEditFlags() {
  MStatus status = MS::kSuccess;
  MArgParser argData = parser();

  // if (argData.isFlagSet(kstrokePointDistFlag)) {
  //     unsigned strokePointDist;
  //     status = argData.getFlagArgument(kstrokePointDistFlag, 0,
  //     strokePointDist);
  // 	if (!status) {
  // 		status.perror("minstrokePointDist flag parsing failed.");
  // 		return status;
  // 	}
  //     fBookmarkContext->setstrokePointDist(strokePointDist);
  // }
  // if (argData.isFlagSet(kStrokeWidthFlag)) {
  //     unsigned strokeWidth;
  //     status = argData.getFlagArgument(kStrokeWidthFlag, 0, strokeWidth);
  // 	if (!status) {
  // 		status.perror("strokeWidth flag parsing failed.");
  // 		return status;
  // 	}
  //     fBookmarkContext->setStrokeWidth(strokeWidth);
  // }

  return status;
}

MStatus BookmarkContextCmd::doQueryFlags() {
  MArgParser argData = parser();
  // if (argData.isFlagSet(kstrokePointDistFlag)) {
  //     setResult((int) fBookmarkContext->getstrokePointDist());
  // }
  // if (argData.isFlagSet(kStrokeWidthFlag)) {
  //     setResult((int) fBookmarkContext->getStrokeWidth());
  // }
  return MS::kSuccess;
}

MPxContext *BookmarkContextCmd::makeObj() {
  // cout << "BookmarkContextCmd::makeObj()" << endl;
  fBookmarkContext = new BookmarkContext();
  return fBookmarkContext;
}

MStatus BookmarkContextCmd::appendSyntax() {
  MSyntax mySyntax = syntax();
  if (MS::kSuccess != mySyntax.addFlag(kStrokePointDistFlag,
                                       kStrokePointDistFlagLong,
                                       MSyntax::kUnsigned)) {
    return MS::kFailure;
  }
  if (MS::kSuccess != mySyntax.addFlag(kStrokeWidthFlag, kStrokeWidthFlagLong,
                                       MSyntax::kUnsigned)) {
    return MS::kFailure;
  }

  return MS::kSuccess;
}

void *BookmarkContextCmd::creator() { return new BookmarkContextCmd; }

} // namespace apps
} // namespace bookmark_squidget