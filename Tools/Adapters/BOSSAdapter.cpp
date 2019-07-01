// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <string>

#include "BOSSAdapter.hpp"
#include "Core/FunctionTreeIntensityWrapper.hpp"
#include "Core/Kinematics.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Physics/ParticleList.hpp"

namespace ComPWA {
namespace Tools {
namespace Adapter {

std::pair<std::shared_ptr<ComPWA::Intensity>, std::shared_ptr<Kinematics>>
BOSS::createHelicityModel(const char *modelXMLFile, int seed,
                          const std::vector<int> &initialState,
                          const std::vector<int> &finalState,
                          const char *particleXMLFile) {
  std::string modelStr;
  if (modelXMLFile) {
    modelStr = std::string(modelXMLFile);
  } else {
    // TODO: error
  }

  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, Physics::defaultParticleList);

  if (particleXMLFile) {
    std::string particleStr(particleXMLFile);
    boost::property_tree::ptree particles;
    boost::property_tree::xml_parser::read_xml(particleXMLFile, particles);
    ReadParticles(partL, particles);
  }

  boost::property_tree::ptree model;
  boost::property_tree::xml_parser::read_xml(modelXMLFile, model);

  ComPWA::Physics::IntensityBuilderXML Builder;
  std::shared_ptr<ComPWA::Intensity> intens;
  std::shared_ptr<Kinematics> kin;
  std::tie(intens, kin) = Builder.createIntensityAndKinematics(model);
  return std::make_pair(intens, kin);
}

} // namespace Adapter
} // namespace Tools
} // namespace ComPWA
