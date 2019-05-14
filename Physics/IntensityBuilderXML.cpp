// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "IntensityBuilderXML.hpp"

#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Physics/CoefficientAmplitudeDecorator.hpp"
#include "Physics/CoherentIntensity.hpp"
#include "Physics/HelicityFormalism/HelicityDecay.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/NormalizationAmplitudeDecorator.hpp"
#include "Physics/NormalizationIntensityDecorator.hpp"
#include "Physics/SequentialAmplitude.hpp"
#include "Physics/StrengthIntensityDecorator.hpp"
#include "Tools/Integration.hpp"

#include "Physics/Dynamics/Flatte.hpp"
#include "Physics/Dynamics/FormFactorDecorator.hpp"
#include "Physics/Dynamics/NonResonant.hpp"
#include "Physics/Dynamics/RelativisticBreitWigner.hpp"
#include "Physics/Dynamics/Voigtian.hpp"

#include <boost/property_tree/ptree.hpp>

#include "ThirdParty/qft++/include/qft++/WignerD.h"

using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

namespace ComPWA {
namespace Physics {

IntensityBuilderXML::IntensityBuilderXML(
    std::shared_ptr<ComPWA::Data::DataSet> phspsample)
    : PhspSample(phspsample) {}

std::tuple<std::shared_ptr<Intensity>, std::shared_ptr<HelicityKinematics>>
IntensityBuilderXML::createIntensityAndKinematics(
    const boost::property_tree::ptree &pt) const {
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, pt);

  auto it = pt.find("HelicityKinematics");
  std::shared_ptr<HelicityKinematics> kin(nullptr);
  if (it != pt.not_found()) {
    kin = createHelicityKinematics(partL, it->second);
  } else {
    throw BadConfig("IntensityBuilderXML::createIntensityAndKinematics(): "
                    " No Kinematics found!");
  }
  it = pt.find("Intensity");
  if (it != pt.not_found())
    return std::make_tuple(createIntensity(partL, kin, it->second), kin);
  else {
    throw BadConfig("IntensityBuilderXML::createIntensityAndKinematics(): "
                    " No Intensity found!");
  }
}

std::shared_ptr<HelicityKinematics>
IntensityBuilderXML::createHelicityKinematics(
    std::shared_ptr<PartList> partL,
    const boost::property_tree::ptree &pt) const {

  auto kininfo = createKinematicsInfo(partL, pt);

  auto phspVal = pt.get_optional<double>("PhspVolume");
  if (phspVal) {
    return std::make_shared<HelicityKinematics>(kininfo, phspVal.get());
  } else {
    return std::make_shared<HelicityKinematics>(kininfo);
  }
}

ParticleStateTransitionKinematicsInfo IntensityBuilderXML::createKinematicsInfo(
    std::shared_ptr<PartList> partL,
    const boost::property_tree::ptree &pt) const {
  auto initialS = pt.get_child("InitialState");
  auto InitialState = std::vector<int>(initialS.size());
  for (auto i : initialS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    unsigned int pos = i.second.get<unsigned int>("<xmlattr>.PositionIndex");
    InitialState.at(pos) = partP.GetId();
  }

  auto finalS = pt.get_child("FinalState");
  auto FinalState = std::vector<int>(finalS.size());
  auto FinalStateEventPositionMapping =
      std::vector<unsigned int>(finalS.size());
  for (auto i : finalS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    unsigned int id = i.second.get<unsigned int>("<xmlattr>.Id");
    unsigned int pos(id);
    boost::optional<unsigned int> opt_pos =
        i.second.get_optional<unsigned int>("<xmlattr>.PositionIndex");
    if (opt_pos)
      pos = opt_pos.get();
    FinalState.at(pos) = partP.GetId();
    FinalStateEventPositionMapping.at(pos) = id;
  }

  if (pt.find("InitialFourMomentum") != pt.not_found()) {
    auto InitialStateP4 =
        createFourMomentum(pt.get_child("InitialFourMomentum"));
    return ParticleStateTransitionKinematicsInfo(
        InitialState, FinalState, partL, InitialStateP4,
        FinalStateEventPositionMapping);
  }

  return ParticleStateTransitionKinematicsInfo(InitialState, FinalState, partL,
                                               FinalStateEventPositionMapping);
}

FourMomentum IntensityBuilderXML::createFourMomentum(
    const boost::property_tree::ptree &pt) const {
  FourMomentum obj;
  double px, py, pz, E;

  auto tmp = pt.get_optional<double>("<xmlattr>.x");
  if (tmp) {
    px = tmp.get();
  } else {
    px = pt.get<double>("x");
  }

  tmp = pt.get_optional<double>("<xmlattr>.y");
  if (tmp) {
    py = tmp.get();
  } else {
    py = pt.get<double>("y");
  }

  tmp = pt.get_optional<double>("<xmlattr>.z");
  if (tmp) {
    pz = tmp.get();
  } else {
    pz = pt.get<double>("z");
  }

  tmp = pt.get_optional<double>("<xmlattr>.E");
  if (tmp) {
    E = tmp.get();
  } else {
    E = pt.get<double>("E");
  }

  obj.setPx(px);
  obj.setPy(py);
  obj.setPz(pz);
  obj.setE(E);
  return obj;
}

std::shared_ptr<Intensity> IntensityBuilderXML::createIntensity(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) const {
  LOG(TRACE) << "loading intensity...";

  std::string IntensityClass(pt.get<std::string>("<xmlattr>.Class"));
  if (IntensityClass == "IncoherentIntensity") {
    return createIncoherentIntensity(partL, kin, pt);
  } else if (IntensityClass == "CoherentIntensity") {
    return createCoherentIntensity(partL, kin, pt);
  } else if (IntensityClass == "StrengthIntensity") {
    return createStrengthIntensity(partL, kin, pt);
  } else if (IntensityClass == "NormalizedIntensity") {
    return createNormalizedIntensity(partL, kin, pt);
  } else {
    throw BadConfig(
        "IntensityBuilderXML::createIntensity() | Found unknown intensity " +
        IntensityClass);
  }
}

std::shared_ptr<Intensity> IntensityBuilderXML::createIncoherentIntensity(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) const {
  std::vector<std::shared_ptr<ComPWA::Intensity>> intensities;
  LOG(TRACE) << "constructing IncoherentIntensity ...";
  auto name = pt.get<std::string>("<xmlattr>.Name");

  for (const auto &x : pt) {
    if (x.first == "Intensity") {
      intensities.push_back(createIntensity(partL, kin, x.second));
    }
  }
  return std::make_shared<ComPWA::Physics::IncoherentIntensity>(name,
                                                                intensities);
}

std::shared_ptr<Intensity> IntensityBuilderXML::createCoherentIntensity(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) const {
  std::vector<std::shared_ptr<ComPWA::Physics::NamedAmplitude>> amps;
  LOG(TRACE) << "constructing CoherentIntensity ...";
  auto name = pt.get<std::string>("<xmlattr>.Name");

  for (const auto &x : pt) {
    if (x.first == "Amplitude") {
      amps.push_back(createAmplitude(partL, kin, x.second));
    }
  }
  return std::make_shared<ComPWA::Physics::CoherentIntensity>(name, amps);
}

std::shared_ptr<Intensity> IntensityBuilderXML::createStrengthIntensity(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) const {
  LOG(TRACE) << "creating StrengthIntensity ...";

  auto name = pt.get<std::string>("<xmlattr>.Name");

  std::shared_ptr<Intensity> UndecoratedIntensity(nullptr);
  std::shared_ptr<FitParameter> Strength(nullptr);

  for (const auto &x : pt) {
    if (x.first == "<xmlattr>")
      continue;
    if (x.first == "Parameter" &&
        x.second.get<std::string>("<xmlattr>.Type") == "Strength") {
      Strength = std::make_shared<FitParameter>(x.second);
    } else if (x.first == "Intensity") {
      UndecoratedIntensity = createIntensity(partL, kin, x.second);
    } else {
      LOG(WARNING) << "IntensityBuilderXML::createStrengthIntensity(): found "
                      "unknown tag "
                   << x.first;
    }
  }

  return std::make_shared<StrengthIntensityDecorator>(
      name, UndecoratedIntensity, Strength);
}

std::shared_ptr<Intensity> IntensityBuilderXML::createNormalizedIntensity(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) const {
  LOG(TRACE) << "creating NormalizedIntensity ...";

  auto name = pt.get<std::string>("<xmlattr>.Name");

  std::shared_ptr<Intensity> UndecoratedIntensity(nullptr);
  std::shared_ptr<Tools::IntegrationStrategy> Integrator(nullptr);
  boost::property_tree::ptree IntegratorPT;

  std::string IntegratorClassName;
  for (const auto &x : pt) {
    if (x.first == "<xmlattr>")
      continue;
    if (x.first == "Intensity") {
      UndecoratedIntensity = createIntensity(partL, kin, x.second);
    } else if (x.first == "IntegrationStrategy") {
      IntegratorPT = x.second;
    } else {
      LOG(WARNING) << "IntensityBuilderXML::createNormalizedIntensity(): found "
                      "unknown tag "
                   << x.first;
    }
  }

  // It is crucial that the IntegrationStrategy is created after the Kinematics
  // is populated with SubSystems from the Intensity
  Integrator = createIntegrationStrategy(partL, kin, IntegratorPT);

  return std::make_shared<NormalizationIntensityDecorator>(
      name, UndecoratedIntensity, Integrator);
}

std::shared_ptr<Tools::IntegrationStrategy>
IntensityBuilderXML::createIntegrationStrategy(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) const {
  LOG(TRACE) << "creating IntegrationStrategy ...";

  if (!pt.empty()) {
    auto ClassName = pt.get<std::string>("<xmlattr>.Class");

    if (ClassName == "MCIntegrationStrategy") {
      if (!PhspSample)
        LOG(FATAL) << "IntensityBuilderXML::createIntegrationStrategy(): phsp "
                      "sample is not set!";
      PhspSample->convertEventsToDataPoints(kin);
      return std::make_shared<ComPWA::Tools::MCIntegrationStrategy>(PhspSample);
    } else {
      LOG(WARNING) << "IntensityBuilderXML::createIntegrationStrategy(): "
                      "IntegrationStrategy type "
                   << ClassName << " unknown!";
    }
  } else {
    LOG(WARNING) << "IntensityBuilderXML::createIntegrationStrategy(): "
                    "IntegrationStrategy tag not specified!";
  }

  LOG(INFO) << "IntensityBuilderXML::createIntegrationStrategy(): creating "
               "default IntegrationStrategy *MCIntegratioStrategy*";
  if (!PhspSample)
    LOG(FATAL) << "IntensityBuilderXML::createIntegrationStrategy(): phsp "
                  "sample is not set!";
  PhspSample->convertEventsToDataPoints(kin);
  return std::make_shared<ComPWA::Tools::MCIntegrationStrategy>(PhspSample);
}

std::shared_ptr<NamedAmplitude> IntensityBuilderXML::createAmplitude(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) const {

  auto ampclass = pt.get<std::string>("<xmlattr>.Class");

  if (ampclass == "HelicityDecay") {
    return createHelicityDecay(
        partL, std::dynamic_pointer_cast<HelicityKinematics>(kin), pt);
  } else if (ampclass == "NormalizedAmplitude") {
    return createNormalizedAmplitude(partL, kin, pt);
  } else if (ampclass == "CoefficientAmplitude") {
    return createCoefficientAmplitude(partL, kin, pt);
  } else if (ampclass == "SequentialAmplitude") {
    return createSequentialAmplitude(partL, kin, pt);
  } else {
    throw BadConfig(
        "IntensityBuilderXML::createAmplitude(): Unknown amplitude " +
        ampclass);
  }
}

std::shared_ptr<NamedAmplitude> IntensityBuilderXML::createNormalizedAmplitude(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) const {
  LOG(TRACE) << "creating NormalizedAmplitude ...";

  auto name = pt.get<std::string>("<xmlattr>.Name");

  std::shared_ptr<NamedAmplitude> UndecoratedAmplitude(nullptr);
  std::shared_ptr<Tools::IntegrationStrategy> Integrator(nullptr);
  boost::property_tree::ptree IntegratorPT;

  std::string IntegratorClassName;
  for (const auto &x : pt) {
    if (x.first == "<xmlattr>")
      continue;
    if (x.first == "Amplitude") {
      UndecoratedAmplitude = createAmplitude(partL, kin, x.second);
    } else if (x.first == "IntegrationStrategy") {
      IntegratorPT = x.second;
    } else {
      LOG(WARNING) << "IntensityBuilderXML::createNormalizedAmplitude(): found "
                      "unknown tag "
                   << x.first;
    }
  }

  // It is crucial that the IntegrationStrategy is created after the Kinematics
  // is populated with SubSystems from the Intensity
  Integrator = createIntegrationStrategy(partL, kin, IntegratorPT);

  return std::make_shared<NormalizationAmplitudeDecorator>(
      name, UndecoratedAmplitude, Integrator);
}

std::shared_ptr<NamedAmplitude> IntensityBuilderXML::createCoefficientAmplitude(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) const {

  LOG(TRACE) << "constructing CoefficientAmplitudeDecorator ...";
  // Name = pt.get<std::string>("<xmlattr>.Name", "empty");

  auto ampname = pt.get<std::string>("<xmlattr>.Name");

  std::shared_ptr<Amplitude> UndecoratedAmplitude(nullptr);
  std::shared_ptr<FitParameter> Magnitude(nullptr);
  std::shared_ptr<FitParameter> Phase(nullptr);
  for (const auto &v : pt) {
    if (v.first == "Parameter") {
      if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude")
        Magnitude = std::make_shared<FitParameter>(v.second);
      if (v.second.get<std::string>("<xmlattr>.Type") == "Phase")
        Phase = std::make_shared<FitParameter>(v.second);
    } else if (v.first == "Amplitude") {
      UndecoratedAmplitude = createAmplitude(partL, kin, v.second);
    }
  }

  if (!Magnitude)
    throw BadParameter(
        "IntensityBuilderXML::createCoefficientAmplitude() | No magnitude "
        "parameter found.");
  if (!Phase)
    throw BadParameter("IntensityBuilderXML::createCoefficientAmplitude() | No "
                       "phase parameter found.");

  return std::make_shared<CoefficientAmplitudeDecorator>(
      ampname, UndecoratedAmplitude, Magnitude, Phase);
}

std::shared_ptr<NamedAmplitude> IntensityBuilderXML::createSequentialAmplitude(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) const {
  LOG(TRACE) << "constructing SequentialAmplitude ...";
  // setName(pt.get<std::string>("<xmlattr>.Name", "empty"));

  auto ampname = pt.get<std::string>("<xmlattr>.Name");

  std::vector<std::shared_ptr<Amplitude>> PartialAmplitudes;
  auto PreFactor = std::complex<double>(1, 0);
  for (const auto &v : pt) {
    if (v.first == "Amplitude") {
      PartialAmplitudes.push_back(createAmplitude(partL, kin, v.second));
    } else if (v.first == "PreFactor") {
      double r = v.second.get<double>("<xmlattr>.Magnitude");
      double p = v.second.get<double>("<xmlattr>.Phase");
      PreFactor = std::polar(r, p);
    } else {
      std::runtime_error(
          "SequentialAmplitude::load() | Cannot create SequentialAmplitude!");
    }
  }
  return std::make_shared<SequentialAmplitude>(ampname, PartialAmplitudes,
                                               PreFactor);
}

std::shared_ptr<NamedAmplitude> IntensityBuilderXML::createHelicityDecay(
    std::shared_ptr<PartList> partL, std::shared_ptr<HelicityKinematics> kin,
    const boost::property_tree::ptree &pt) const {

  LOG(TRACE) << "HelicityDecay::load() |";
  unsigned int SubSystemIndex(kin->addSubSystem(SubSystem(pt)));
  // SubSys = kin->subSystem(SubSystemIndex);
  unsigned int DataPosition = 3 * SubSystemIndex;

  auto ampname = pt.get<std::string>("<xmlattr>.Name");
  // Name = pt.get<std::string>("<xmlattr>.Name", "empty");

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  auto partItr = partL->find(name);
  if (partItr == partL->end())
    throw std::runtime_error("HelicityDecay::load | Particle " + name +
                             " not found in list!");
  ComPWA::Spin J = partItr->second.GetSpinQuantumNumber("Spin");
  ComPWA::Spin mu(pt.get<double>("DecayParticle.<xmlattr>.Helicity"));
  // if the node OrbitalAngularMomentum does not exist, set it to spin J as
  // default value
  ComPWA::Spin orbitL(pt.get<double>(
      "DecayParticle.<xmlattr>.OrbitalAngularMomentum", (double)J));

  double PreFactor(1.0);
  const auto &canoSum = pt.get_child_optional("CanonicalSum");
  if (canoSum) {
    const auto &sumTree = canoSum.get();
    orbitL = sumTree.get<double>("<xmlattr>.L");
    double coef = sqrt((2 * (double)orbitL + 1) / (2 * (double)J + 1));
    for (const auto &cg : sumTree.get_child("")) {
      if (cg.first != "ClebschGordan")
        continue;
      double j1 = cg.second.get<double>("<xmlattr>.j1");
      double m1 = cg.second.get<double>("<xmlattr>.m1");
      double j2 = cg.second.get<double>("<xmlattr>.j2");
      double m2 = cg.second.get<double>("<xmlattr>.m2");
      double J = cg.second.get<double>("<xmlattr>.J");
      double M = cg.second.get<double>("<xmlattr>.M");
      coef *= ComPWA::QFT::Clebsch(j1, m1, j2, m2, J, M);
    }
    PreFactor *= coef;
  }

  // Read name and helicities from decay products
  auto decayProducts = pt.get_child("DecayProducts");
  if (decayProducts.size() != 2)
    throw boost::property_tree::ptree_error(
        "IntensityBuilderXML::createHelicityDecay() | Expect exactly two decay "
        "products (" +
        std::to_string(decayProducts.size()) + " given)!");

  std::pair<std::string, std::string> DecayProducts;
  std::pair<ComPWA::Spin, ComPWA::Spin> DecayHelicities;

  auto p = decayProducts.begin();
  DecayProducts.first = p->second.get<std::string>("<xmlattr>.Name");
  DecayHelicities.first =
      ComPWA::Spin(p->second.get<double>("<xmlattr>.Helicity"));
  ++p;
  DecayProducts.second = p->second.get<std::string>("<xmlattr>.Name");
  DecayHelicities.second =
      ComPWA::Spin(p->second.get<double>("<xmlattr>.Helicity"));

  auto partProp = partL->find(name)->second;
  std::string decayType = partProp.GetDecayType();

  std::shared_ptr<ComPWA::Physics::Dynamics::AbstractDynamicalFunction>
      DynamicFunction(nullptr);

  if (decayType == "stable") {
    throw std::runtime_error(
        "IntensityBuilderXML::createHelicityDecay() | Stable particle is "
        "given as mother particle of a decay. Makes no sense!");
  } else if (decayType == "relativisticBreitWigner") {
    using ComPWA::Physics::Dynamics::RelativisticBreitWigner;
    auto RBW = new RelativisticBreitWigner(name, DecayProducts, partL);
    RBW->SetOrbitalAngularMomentum(orbitL);
    DynamicFunction = std::shared_ptr<RelativisticBreitWigner>(RBW);
  } else if (decayType == "flatte") {
    using ComPWA::Physics::Dynamics::Flatte;
    auto F = new Flatte(name, DecayProducts, partL);
    F->SetOrbitalAngularMomentum(orbitL);
    DynamicFunction = std::shared_ptr<Flatte>(F);
  } else if (decayType == "voigt") {
    DynamicFunction =
        std::make_shared<Dynamics::Voigtian>(name, DecayProducts, partL);
  } else if (decayType == "virtual" || decayType == "nonResonant") {
    DynamicFunction = std::make_shared<Dynamics::NonResonant>(name);
  } else {
    throw std::runtime_error("HelicityDecay::Factory() | Unknown decay type " +
                             decayType + "!");
  }

  // set production formfactor
  std::string daug1Name = DecayProducts.first;
  std::string daug2Name = DecayProducts.second;
  auto parMass1 = std::make_shared<FitParameter>(
      partL->find(daug1Name)->second.GetMassPar());
  auto parMass2 = std::make_shared<FitParameter>(
      partL->find(daug2Name)->second.GetMassPar());
  auto decayInfo = partProp.GetDecayInfo();
  int ffType = 0;
  std::shared_ptr<ComPWA::FitParameter> parRadius;
  for (const auto &node : decayInfo.get_child("")) {
    if (node.first == "FormFactor") {
      ffType = node.second.get<int>("<xmlattr>.Type");
    } else if (node.first == "Parameter") {
      std::string parType = node.second.get<std::string>("<xmlattr>.Type");
      if (parType == "MesonRadius") {
        parRadius = std::make_shared<ComPWA::FitParameter>(node.second);
      }
    }
  }

  std::shared_ptr<HelicityFormalism::HelicityDecay> HeliDecay;

  if (ffType == 0 || ((unsigned int)orbitL == 0)) {
    HeliDecay = std::make_shared<HelicityFormalism::HelicityDecay>(
        ampname,
        std::make_shared<HelicityFormalism::AmpWignerD>(
            J, mu, DecayHelicities.first - DecayHelicities.second),
        DynamicFunction, DataPosition, PreFactor);
  } else {
    if (parRadius == nullptr) {
      throw std::runtime_error(
          "IntensityBuilderXML::createHelicityDecay() | no MesonRadius is "
          "given! It is needed to calculate the formfactor!");
    }

    std::shared_ptr<ComPWA::Physics::Dynamics::AbstractDynamicalFunction>
        DyFuncWithProductionFF =
            std::make_shared<ComPWA::Physics::Dynamics ::FormFactorDecorator>(
                name, DynamicFunction, parMass1, parMass2, parRadius, orbitL,
                (ComPWA::Physics::Dynamics::FormFactorType)ffType);

    HeliDecay = std::make_shared<HelicityFormalism::HelicityDecay>(
        ampname,
        std::make_shared<HelicityFormalism::AmpWignerD>(
            J, mu, DecayHelicities.first - DecayHelicities.second),
        DynamicFunction, DataPosition, PreFactor);
  }

  return HeliDecay;
}

} // namespace Physics
} // namespace ComPWA
