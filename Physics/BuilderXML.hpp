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
#include "Physics/SubSystem.hpp"
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
  IntensityBuilderXML(ParticleList ParticleList, Kinematics &Kinematics,
                      const boost::property_tree::ptree &ModelTree,
                      const EventCollection &TruePhspSample = {},
                      const EventCollection &RecoPhspSample = {});

  ComPWA::FunctionTree::FunctionTreeIntensity createIntensity();

  std::vector<ComPWA::Tools::IntensityComponent> createIntensityComponents(
      std::vector<std::vector<std::string>> ComponentList = {});

  std::map<std::string, std::string> getAllComponentNames() const;

private:
  struct DataContainer {
    ComPWA::FunctionTree::ParameterList Data;
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>> Weights;
    double WeightSum = 0.0;
  };

  std::shared_ptr<ComPWA::FunctionTree::TreeNode>
  createIntensityFT(const boost::property_tree::ptree &pt,
                    const ComPWA::FunctionTree::ParameterList &DataSample);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> createIncoherentIntensityFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> createIncoherentIntensityFT(
      std::vector<std::shared_ptr<ComPWA::FunctionTree::TreeNode>> Intensities);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> createCoherentIntensityFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> createCoherentIntensityFT(
      std::vector<std::shared_ptr<ComPWA::FunctionTree::TreeNode>> Amplitudes);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> createStrengthIntensityFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> createNormalizedIntensityFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode>
  normalizeIntensityFT(const boost::property_tree::ptree &UnnormalizedPT,
                       const ComPWA::FunctionTree::ParameterList &DataSample,
                       std::string IntegratorClassNa);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> createIntegrationStrategyFT(
      std::shared_ptr<ComPWA::FunctionTree::TreeNode> UnnormalizedIntensity,
      std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
          PhspWeights,
      double PhspWeightSum, std::string IntegratorClassName);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode>
  createAmplitudeFT(const boost::property_tree::ptree &pt,
                    const ComPWA::FunctionTree::ParameterList &DataSample);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> createNormalizedAmplitudeFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> createCoefficientAmplitudeFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> createSequentialAmplitudeFT(
      const boost::property_tree::ptree &pt,
      const ComPWA::FunctionTree::ParameterList &DataSample);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode>
  createHelicityDecayFT(const boost::property_tree::ptree &pt,
                        const ComPWA::FunctionTree::ParameterList &DataSample);

  void updateDataContainerWeights(DataContainer &DataCon,
                                  const EventCollection &DataSample);
  void updateDataContainerState();
  void updateDataContainerContent();

  void
  addFunctionTreeComponent(std::string Name, std::string Type,
                           std::shared_ptr<ComPWA::FunctionTree::TreeNode> FT);

  bool ComponentRegisteringEnabled = true;

  std::map<
      std::string,
      std::pair<std::string, std::shared_ptr<ComPWA::FunctionTree::TreeNode>>>
      UniqueComponentFTMapping;

  ParticleList PartList;
  Kinematics &Kinematic;
  boost::property_tree::ptree ModelTree;
  const EventCollection &TruePhspSample;
  const EventCollection &RecoPhspSample;

  ComPWA::FunctionTree::ParameterList Parameters;
  DataContainer Data;
  DataContainer PhspData;
  DataContainer PhspRecoData;
};

/// Create HelicityKinematics object from an XML file that contains both a
/// kinematics section and a particle section
HelicityFormalism::HelicityKinematics
createHelicityKinematics(const std::string XmlFile);

/// Create HelicityKinematics object from an XML file with a kinematics section
/// and provide a particle list separately.
HelicityFormalism::HelicityKinematics
createHelicityKinematics(const ComPWA::ParticleList &PartList,
                         const std::string XmlFile);

HelicityFormalism::HelicityKinematics
createHelicityKinematics(const ComPWA::ParticleList &PartList,
                         const boost::property_tree::ptree &pt);

struct TwoBodyDecayInfo {
  SubSystem SubSys;
  std::pair<std::string, std::string> Names;
  std::pair<double, double> Helicities;
};

TwoBodyDecayInfo extractDecayInfo(const boost::property_tree::ptree &pt);

} // namespace Physics
} // namespace ComPWA

#endif
