#include "bookmark_squidget/canvas/attributeConfigManager.hpp"

namespace bookmark_squidget {
namespace canvas {

AttributeConfigManager::AttributeConfigManager() { interpolationType = LINEAR; }

AttributeConfigManager::~AttributeConfigManager() {}

void AttributeConfigManager::addConfig(
    std::shared_ptr<AttrConfigMap> configMap) {
  /**
   * @brief Adds new keyed configuration to the attributeConfigSteps vector at
   * time t.
   */
  AttributeConfigStep step;
  step.attrValueMap = configMap;
  cout << "configMap->size(): " << configMap->size() << endl;
  cout << "step.attrValueMap->size(): " << step.attrValueMap->size() << endl;
  keyedConfigs.push_back(step);
}

void AttributeConfigManager::removeAttributeConfigStep(int index) {
  /**
   * @brief Deletes keyed value at index.
   */
  keyedConfigs.erase(keyedConfigs.begin() + index);
}

AttrConfigMap AttributeConfigManager::interpolate(double t) {
  /**
   * @brief Find interpolated attribute values for a given time t.
   */
  // cout << "AttributeConfigManager::interpolate()" << endl;
  AttributeConfigStep a = keyedConfigs[0];
  AttributeConfigStep b = keyedConfigs[1];

  // cout << keyedConfigs.size() << endl;
  // cout << a.stepValue << " " << b.stepValue << endl;
  int index = 1;
  while (t > b.stepValue) {
    a = b;
    b = keyedConfigs[++index];
  }

  AttrConfigMap interpolatedMap;
  std::ostringstream ss;
  for (auto &pair : *(a.attrValueMap)) {
    double aVal = pair.second;
    double bVal = b.attrValueMap->at(pair.first);
    double interpolatedVal =
        aVal + (bVal - aVal) * (t - a.stepValue) / (b.stepValue - a.stepValue);
    interpolatedMap[pair.first] = interpolatedVal;
  }

  // cout << "command: " << ss.str() << endl;
  MGlobal::executeCommand(ss.str().c_str());

  return interpolatedMap;
}

void AttributeConfigManager::sortAttributeConfigSteps() {
  /**
   * @brief Sort the attributeConfigSteps vector by stepValue.
   */
  // cout << "sortAttributeConfigSteps() is disabled" << endl;
  // sort(begin(keyedConfigs), end(keyedConfigs),
  //      [](AttributeConfigStep a, AttributeConfigStep b) {
  //        return a.stepValue < b.stepValue;
  //      });
}

std::vector<AttributeConfigStep>
AttributeConfigManager::getChangedConfigs(bool isSingleStroke) {
  std::vector<AttributeConfigStep> configs;

  // Initialize new, empty config vector with correct step values.
  for (size_t i = 0; i < keyedConfigs.size(); ++i) {
    AttributeConfigStep newConfig;
    newConfig.stepValue = keyedConfigs.at(i).stepValue;
    newConfig.attrValueMap = std::make_shared<AttrConfigMap>();

    configs.push_back(newConfig);
  }

  std::set<std::string> changingPlugNames;

  // For all the keyed configs, find the plugs that change
  for (size_t i = 0; i < keyedConfigs.size() - 1; ++i) {
    auto currConfigMap = keyedConfigs.at(i + 0).attrValueMap;
    auto nextConfigMap = keyedConfigs.at(i + 1).attrValueMap;

    // Check if key in currConfigMap is in nextConfigMap
    for (auto const &configPair : *currConfigMap) {
      std::string currPlugName = configPair.first;
      double currPlugValue = configPair.second;

      auto it = nextConfigMap->find(currPlugName);
      if (it == nextConfigMap->end()) {
        break;
      }

      double nextPlugValue = nextConfigMap->at(currPlugName);
      bool plugChanged = currPlugValue != nextPlugValue;
      if (plugChanged || isSingleStroke) {
        changingPlugNames.emplace(currPlugName);
      }
    }
  }

  for (size_t i = 0; i < keyedConfigs.size(); ++i) {
    auto currConfigMap = keyedConfigs.at(i).attrValueMap;
    AttributeConfigStep &currConfigStep = configs.at(i);

    for (auto const plugName : changingPlugNames) {
      double val = currConfigMap->at(plugName);

      if (currConfigMap->find(plugName) != currConfigMap->end()) {
        currConfigStep.attrValueMap->emplace(plugName, val);
      }
    }
  }

  // cout << "Ending configs.size(): " << configs.size() << endl;
  // for (size_t i = 0; i < keyedConfigs.size(); ++i) {
  //   AttributeConfigStep &currConfigStep = configs.at(i);
  //   for (auto const &pair : *(currConfigStep.attrValueMap)) {
  // cout << "currConfigStep: " << i << " " << pair.first << " " << pair.second
  // << endl;
  //   }
  // }

  // for (size_t i = 0; i < configs.size(); ++i) {
  //   AttributeConfigStep &currConfigStep = configs.at(i);
  //   // cout << "currConfigStep: " << i << " " << currConfigStep.stepValue <<
  //   endl;
  // }
  return configs;
}

} // namespace canvas
} // namespace bookmark_squidget