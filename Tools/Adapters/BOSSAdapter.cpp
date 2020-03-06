// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <string>

#include "BOSSAdapter.hpp"
#include "Core/Kinematics.hpp"
#include "Physics/BuilderXML.hpp"

#include "boost/property_tree/xml_parser.hpp"

namespace ComPWA {
namespace Tools {
namespace Adapter {

std::pair<FunctionTree::FunctionTreeIntensity,
          Physics::HelicityFormalism::HelicityKinematics>
BOSS::createHelicityModel(const char *modelXMLFile, int seed,
                          const std::vector<pid> &initialState,
                          const std::vector<pid> &finalState,
                          const char *particleXMLFile) {
  std::string modelStr;
  if (modelXMLFile == nullptr) {
    modelStr = std::string(modelXMLFile);
  } else {
    // TODO: error
  }

  ComPWA::ParticleList Particles;

  if (particleXMLFile) {
    std::string particleStr(particleXMLFile);
    boost::property_tree::ptree ParticlesPT;
    boost::property_tree::xml_parser::read_xml(particleXMLFile, ParticlesPT);
    ComPWA::insertParticles(Particles, ParticlesPT);
  }

  boost::property_tree::ptree model;
  boost::property_tree::xml_parser::read_xml(modelXMLFile, model);

  auto kin(ComPWA::Physics::createHelicityKinematics(Particles, model));

  ComPWA::Physics::IntensityBuilderXML Builder(Particles, kin, model);

  return std::make_pair(Builder.createIntensity(), std::move(kin));
}

} // namespace Adapter
} // namespace Tools
} // namespace ComPWA
