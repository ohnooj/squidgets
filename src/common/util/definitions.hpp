#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

#include <map>
#include <string>

#include <Eigen/Dense>

#include <maya/MPlugArray.h>
#include <maya/MPointArray.h>
#include <maya/MString.h>

namespace common {
namespace util {

typedef std::string MapKey;
typedef std::map<MapKey, MPlugArray> AttrMap;
typedef std::map<MapKey, Eigen::VectorXd> AttrValueMap;
typedef MPointArray PenPoints;

} // namespace util
} // namespace common
#endif