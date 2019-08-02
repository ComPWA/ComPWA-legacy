// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_INTENSITYBUILDER_HPP_
#define COMPWA_PHYSICS_INTENSITYBUILDER_HPP_

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

struct IntensityBuilderState {
  ComPWA::FunctionTree::ParameterList Parameters;
  ComPWA::FunctionTree::ParameterList Data;
  ComPWA::FunctionTree::ParameterList PhspData;
  ComPWA::FunctionTree::ParameterList ActiveData;
  bool IsDataActive = true;
};

class IntensityBuilderXML {
public:
  IntensityBuilderXML(const std::vector<Event> phspsample = {});

  std::pair<ComPWA::FunctionTree::FunctionTreeIntensity,
            HelicityFormalism::HelicityKinematics>
  createIntensityAndKinematics(const boost::property_tree::ptree &pt);

  HelicityFormalism::HelicityKinematics
  createHelicityKinematics(std::shared_ptr<PartList> partL,
                           const boost::property_tree::ptree &pt);

  ComPWA::FunctionTree::FunctionTreeIntensity
  createIntensity(std::shared_ptr<PartList> partL, Kinematics &kin,
                  const boost::property_tree::ptree &pt);

private:
  ParticleStateTransitionKinematicsInfo
  createKinematicsInfo(std::shared_ptr<PartList> partL,
                       const boost::property_tree::ptree &pt) const;

  FourMomentum createFourMomentum(const boost::property_tree::ptree &pt) const;

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createIntensityFT(std::shared_ptr<PartList> partL, Kinematics &kin,
                    const boost::property_tree::ptree &pt, std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createIncoherentIntensityFT(std::shared_ptr<PartList> partL, Kinematics &kin,
                              const boost::property_tree::ptree &pt,
                              std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createCoherentIntensityFT(std::shared_ptr<PartList> partL, Kinematics &kin,
                            const boost::property_tree::ptree &pt,
                            std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createStrengthIntensityFT(std::shared_ptr<PartList> partL, Kinematics &kin,
                            const boost::property_tree::ptree &pt,
                            std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createNormalizedIntensityFT(std::shared_ptr<PartList> partL, Kinematics &kin,
                              const boost::property_tree::ptree &pt,
                              std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  normalizeIntensityFT(std::shared_ptr<PartList> partL, Kinematics &kin,
                       const boost::property_tree::ptree &UnnormalizedPT,
                       std::string IntegratorClassName);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createIntegrationStrategyFT(
      std::shared_ptr<ComPWA::FunctionTree::FunctionTree> UnnormalizedIntensity,
      Kinematics &kin, std::string IntegratorClassName);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createAmplitudeFT(std::shared_ptr<PartList> partL, Kinematics &kin,
                    const boost::property_tree::ptree &pt, std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createNormalizedAmplitudeFT(std::shared_ptr<PartList> partL, Kinematics &kin,
                              const boost::property_tree::ptree &pt);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createCoefficientAmplitudeFT(std::shared_ptr<PartList> partL, Kinematics &kin,
                               const boost::property_tree::ptree &pt,
                               std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createSequentialAmplitudeFT(std::shared_ptr<PartList> partL, Kinematics &kin,
                              const boost::property_tree::ptree &pt,
                              std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createHelicityDecayFT(std::shared_ptr<PartList> partL,
                        HelicityFormalism::HelicityKinematics &kin,
                        const boost::property_tree::ptree &pt,
                        std::string suffix);

  void updateDataContainerState(const Kinematics &kin);

  std::vector<ComPWA::Event> PhspSample;
  std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>> PhspWeights;

  IntensityBuilderState CurrentIntensityState;
};

} // namespace Physics
} // namespace ComPWA

#endif
