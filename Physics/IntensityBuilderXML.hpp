// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_INTENSITYBUILDER_HPP_
#define COMPWA_PHYSICS_INTENSITYBUILDER_HPP_

#include <memory>
#include <tuple>
#include <vector>

#include "Core/Function.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"

#include <boost/property_tree/ptree_fwd.hpp>

namespace ComPWA {
class Kinematics;
class OldIntensity;
class FunctionTreeIntensityWrapper;

namespace Data {
class DataSet;
}

namespace Tools {
class IntegrationStrategy;
}

namespace Physics {
class NamedAmplitude;

namespace HelicityFormalism {
class HelicityKinematics;
}

class IntensityBuilderXML {
public:
  IntensityBuilderXML(std::shared_ptr<ComPWA::Data::DataSet> phspsample = {});

  std::tuple<std::shared_ptr<FunctionTreeIntensityWrapper>,
             std::shared_ptr<HelicityFormalism::HelicityKinematics>>
  createIntensityAndKinematics(const boost::property_tree::ptree &pt) const;

  std::shared_ptr<HelicityFormalism::HelicityKinematics>
  createHelicityKinematics(std::shared_ptr<PartList> partL,
                           const boost::property_tree::ptree &pt) const;
  std::shared_ptr<ComPWA::FunctionTreeIntensityWrapper>
  createIntensity(std::shared_ptr<PartList> partL,
                  std::shared_ptr<Kinematics> kin,
                  const boost::property_tree::ptree &pt) const;

  std::shared_ptr<ComPWA::OldIntensity>
  createOldIntensity(std::shared_ptr<PartList> partL,
                     std::shared_ptr<Kinematics> kin,
                     const boost::property_tree::ptree &pt) const;

  std::shared_ptr<NamedAmplitude>
  createAmplitude(std::shared_ptr<PartList> partL,
                  std::shared_ptr<Kinematics> kin,
                  const boost::property_tree::ptree &pt) const;

private:
  ParticleStateTransitionKinematicsInfo
  createKinematicsInfo(std::shared_ptr<PartList> partL,
                       const boost::property_tree::ptree &pt) const;

  FourMomentum createFourMomentum(const boost::property_tree::ptree &pt) const;

  std::shared_ptr<OldIntensity>
  createIncoherentIntensity(std::shared_ptr<PartList> partL,
                            std::shared_ptr<Kinematics> kin,
                            const boost::property_tree::ptree &pt) const;
  std::shared_ptr<OldIntensity>
  createCoherentIntensity(std::shared_ptr<PartList> partL,
                          std::shared_ptr<Kinematics> kin,
                          const boost::property_tree::ptree &pt) const;

  std::shared_ptr<OldIntensity>
  createStrengthIntensity(std::shared_ptr<PartList> partL,
                          std::shared_ptr<Kinematics> kin,
                          const boost::property_tree::ptree &pt) const;

  std::shared_ptr<OldIntensity>
  createNormalizedIntensity(std::shared_ptr<PartList> partL,
                            std::shared_ptr<Kinematics> kin,
                            const boost::property_tree::ptree &pt) const;

  std::shared_ptr<Tools::IntegrationStrategy>
  createIntegrationStrategy(std::shared_ptr<PartList> partL,
                            std::shared_ptr<Kinematics> kin,
                            const boost::property_tree::ptree &pt) const;

  std::shared_ptr<NamedAmplitude>
  createNormalizedAmplitude(std::shared_ptr<PartList> partL,
                            std::shared_ptr<Kinematics> kin,
                            const boost::property_tree::ptree &pt) const;

  std::shared_ptr<NamedAmplitude>
  createCoefficientAmplitude(std::shared_ptr<PartList> partL,
                             std::shared_ptr<Kinematics> kin,
                             const boost::property_tree::ptree &pt) const;

  std::shared_ptr<NamedAmplitude>
  createSequentialAmplitude(std::shared_ptr<PartList> partL,
                            std::shared_ptr<Kinematics> kin,
                            const boost::property_tree::ptree &pt) const;

  std::shared_ptr<NamedAmplitude> createHelicityDecay(
      std::shared_ptr<PartList> partL,
      std::shared_ptr<HelicityFormalism::HelicityKinematics> kin,
      const boost::property_tree::ptree &pt) const;

  std::shared_ptr<ComPWA::Data::DataSet> PhspSample;
};

} // namespace Physics
} // namespace ComPWA

#endif
