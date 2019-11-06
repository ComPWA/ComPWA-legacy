// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_BUILDERXML_HPP_
#define COMPWA_PHYSICS_BUILDERXML_HPP_

#include "Core/Function.hpp"
#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/FunctionTree/Value.hpp"
#include "Data/DataSet.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"
#include "Tools/FitFractions.hpp"

#include <boost/property_tree/ptree_fwd.hpp>

#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

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
  IntensityBuilderXML(ParticleList PartList_, Kinematics &Kin,
                      const boost::property_tree::ptree &ModelTree_,
                      const std::vector<Event> &TruePhspSample_ = {},
                      const std::vector<Event> &RecoPhspSample_ = {});

  ComPWA::FunctionTree::FunctionTreeIntensity createIntensity();

  std::vector<ComPWA::Tools::IntensityComponent> createIntensityComponents(
      std::vector<std::vector<std::string>> ComponentList = {});

private:
  struct DataContainer {
    ComPWA::FunctionTree::ParameterList Data;
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>> Weights;
    double WeightSum = 0.0;
  };

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createIntensityFT(const boost::property_tree::ptree &pt,
                    const ComPWA::FunctionTree::ParameterList &DataSample,
                    std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createIncoherentIntensityFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample,
      std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createIncoherentIntensityFT(
      std::string Name,
      std::vector<std::shared_ptr<ComPWA::FunctionTree::FunctionTree>>
          Intensities,
      std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> createCoherentIntensityFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample,
      std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> createCoherentIntensityFT(
      std::string Name,
      std::vector<std::shared_ptr<ComPWA::FunctionTree::FunctionTree>>
          Amplitudes,
      std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> createStrengthIntensityFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample,
      std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createNormalizedIntensityFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample,
      std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  normalizeIntensityFT(const boost::property_tree::ptree &UnnormalizedPT,
                       const ComPWA::FunctionTree::ParameterList &DataSample,
                       std::string IntegratorClassName, std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createIntegrationStrategyFT(
      std::shared_ptr<ComPWA::FunctionTree::FunctionTree> UnnormalizedIntensity,
      std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
          PhspWeights,
      double PhspWeightSum, std::string IntegratorClassName);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createAmplitudeFT(const boost::property_tree::ptree &pt,
                    const ComPWA::FunctionTree::ParameterList &DataSample,
                    std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createNormalizedAmplitudeFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample,
      std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createCoefficientAmplitudeFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample,
      std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createSequentialAmplitudeFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample,
      std::string suffix);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createHelicityDecayFT(const boost::property_tree::ptree &pt,
                        const ComPWA::FunctionTree::ParameterList &DataSample,
                        std::string suffix);

  void updateDataContainerWeights(DataContainer &DataCon,
                                  const std::vector<ComPWA::Event> &DataSample);
  void updateDataContainerState();
  void updateDataContainerContent();

  void addFunctionTreeComponent(
      std::string Name, std::string Type,
      std::shared_ptr<ComPWA::FunctionTree::FunctionTree> FT);

  std::map<std::string,
           std::pair<std::string,
                     std::shared_ptr<ComPWA::FunctionTree::FunctionTree>>>
      UniqueComponentFTMapping;

  ParticleList PartList;
  Kinematics &Kinematic;
  boost::property_tree::ptree ModelTree;
  const std::vector<ComPWA::Event> &TruePhspSample;
  const std::vector<ComPWA::Event> &RecoPhspSample;

  ComPWA::FunctionTree::ParameterList Parameters;
  DataContainer Data;
  DataContainer PhspData;
  DataContainer PhspRecoData;
};

HelicityFormalism::HelicityKinematics
createHelicityKinematics(const ComPWA::ParticleList &PartList,
                         const boost::property_tree::ptree &pt);

} // namespace Physics
} // namespace ComPWA

#endif
