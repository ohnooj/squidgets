#ifndef ATTRIBUTE_CONFIG_MANAGER_HPP
#define ATTRIBUTE_CONFIG_MANAGER_HPP

#include <map>
#include <memory>
#include <set>
#include <vector>

#include <maya/MGlobal.h>
#include <maya/MIntArray.h>

#include "common/util/definitions.hpp"

namespace bookmark_squidget {
namespace canvas {
namespace {
using AttrConfigMap = std::map<std::string, double>;
}

enum InterpolationType { LINEAR, CUBIC };

struct AttributeConfigStep {
  double stepValue;
  std::shared_ptr<AttrConfigMap> attrValueMap;
};

class AttributeConfigManager {
public:
  AttributeConfigManager();
  ~AttributeConfigManager();

  void addConfig(std::shared_ptr<AttrConfigMap> configMap);
  void removeAttributeConfigStep(int index);

  AttrConfigMap interpolate(double t);
  std::vector<AttributeConfigStep> getChangedConfigs(bool isSingleStroke);

private:
  InterpolationType interpolationType;
  std::vector<AttributeConfigStep> keyedConfigs; // size must be >= 2

  void sortAttributeConfigSteps();
};
} // namespace canvas
} // namespace bookmark_squidget
#endif
