// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_BUILDERXML_HPP_
#define COMPWA_PHYSICS_BUILDERXML_HPP_

#include <memory>
#include <tuple>
#include <vector>

#include "Core/Function.hpp"
#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/FunctionTree/Value.hpp"
#include "Data/DataSet.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"

#include <boost/property_tree/ptree_fwd.hpp>

namespace ComPWA {
class Kinematics;

namespace Tools {
class IntegrationStrategy;
}

namespace Physics {
namespace HelicityFormalism {
class HelicityKinematics;
}

class IntensityBuilderXML {
public:
  IntensityBuilderXML(std::shared_ptr<PartList> ParticlList_, Kinematics &Kin,
                      const boost::property_tree::ptree &ModelTree_,
                      std::vector<Event> PhspSample_ = {});

  ComPWA::FunctionTree::FunctionTreeIntensity createIntensity();

private:
  struct IntensityBuilderState {
    ComPWA::FunctionTree::ParameterList Parameters;
    ComPWA::FunctionTree::ParameterList Data;
    ComPWA::FunctionTree::ParameterList PhspData;
    ComPWA::FunctionTree::ParameterList ActiveData;
    bool IsDataActive = true;
  };

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createIntensityFT(const boost::property_tree::ptree &pt, std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createIncoherentIntensityFT(const boost::property_tree::ptree &pt,
                              std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createCoherentIntensityFT(const boost::property_tree::ptree &pt,
                            std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createStrengthIntensityFT(const boost::property_tree::ptree &pt,
                            std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createNormalizedIntensityFT(const boost::property_tree::ptree &pt,
                              std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  normalizeIntensityFT(const boost::property_tree::ptree &UnnormalizedPT,
                       std::string IntegratorClassName);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createIntegrationStrategyFT(
      std::shared_ptr<ComPWA::FunctionTree::FunctionTree> UnnormalizedIntensity,
      std::string IntegratorClassName);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createAmplitudeFT(const boost::property_tree::ptree &pt, std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createNormalizedAmplitudeFT(const boost::property_tree::ptree &pt);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createCoefficientAmplitudeFT(const boost::property_tree::ptree &pt,
                               std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createSequentialAmplitudeFT(const boost::property_tree::ptree &pt,
                              std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createHelicityDecayFT(const boost::property_tree::ptree &pt,
                        std::string suffix);

  void updateDataContainerState();

  std::shared_ptr<PartList> ParticleList;
  Kinematics &Kinematic;
  boost::property_tree::ptree ModelTree;
  std::vector<ComPWA::Event> PhspSample;
  std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>> PhspWeights;

  IntensityBuilderState CurrentIntensityState;
};

HelicityFormalism::HelicityKinematics
createHelicityKinematics(std::shared_ptr<PartList> partL,
                         const boost::property_tree::ptree &pt);

ParticleStateTransitionKinematicsInfo
createKinematicsInfo(std::shared_ptr<PartList> partL,
                     const boost::property_tree::ptree &pt);

FourMomentum createFourMomentum(const boost::property_tree::ptree &pt);

} // namespace Physics
} // namespace ComPWA

#endif
