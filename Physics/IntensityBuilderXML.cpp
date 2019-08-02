// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "IntensityBuilderXML.hpp"

#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Tools/Integration.hpp"

#include "Physics/Dynamics/Flatte.hpp"
#include "Physics/Dynamics/FormFactor.hpp"
#include "Physics/Dynamics/NonResonant.hpp"
#include "Physics/Dynamics/RelativisticBreitWigner.hpp"
#include "Physics/Dynamics/Voigtian.hpp"
#include "Physics/HelicityFormalism/WignerD.hpp"

#include <boost/property_tree/ptree.hpp>

#include "ThirdParty/qft++/include/qft++/WignerD.h"

using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

namespace ComPWA {
namespace Physics {

using ComPWA::FunctionTree::FitParameter;
using ComPWA::FunctionTree::FunctionTreeIntensity;

IntensityBuilderXML::IntensityBuilderXML(std::vector<Event> phspsample)
    : PhspSample(phspsample) {}

std::pair<FunctionTreeIntensity, HelicityKinematics>
IntensityBuilderXML::createIntensityAndKinematics(
    const boost::property_tree::ptree &pt) {
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, pt);

  auto it = pt.find("HelicityKinematics");
  if (it != pt.not_found()) {
    HelicityKinematics kin(createHelicityKinematics(partL, it->second));

    it = pt.find("Intensity");
    if (it != pt.not_found()) {
      auto intens = createIntensity(partL, kin, it->second);
      return std::make_pair(intens, std::move(kin));
    } else {
      throw BadConfig("IntensityBuilderXML::createIntensityAndKinematics(): "
                      " No Intensity found!");
    }
  } else {
    throw BadConfig("IntensityBuilderXML::createIntensityAndKinematics(): "
                    " No Kinematics found!");
  }
}

HelicityKinematics IntensityBuilderXML::createHelicityKinematics(
    std::shared_ptr<PartList> partL, const boost::property_tree::ptree &pt) {

  auto kininfo = createKinematicsInfo(partL, pt);

  auto phspVal = pt.get_optional<double>("PhspVolume");
  if (phspVal) {
    return HelicityKinematics(kininfo, phspVal.get());
  } else {
    return HelicityKinematics(kininfo);
  }
}

ComPWA::FunctionTree::FunctionTreeIntensity
IntensityBuilderXML::createIntensity(std::shared_ptr<PartList> partL,
                                     Kinematics &kin,
                                     const boost::property_tree::ptree &pt) {
  LOG(TRACE) << "loading intensity...";

  CurrentIntensityState = IntensityBuilderState(); // BlankState

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> FT =
      createIntensityFT(partL, kin, pt, "");

  CurrentIntensityState.Data = CurrentIntensityState.ActiveData;

  return ComPWA::FunctionTree::FunctionTreeIntensity(
      FT, CurrentIntensityState.Parameters, CurrentIntensityState.Data);
}

ParticleStateTransitionKinematicsInfo IntensityBuilderXML::createKinematicsInfo(
    std::shared_ptr<PartList> partL,
    const boost::property_tree::ptree &pt) const {
  auto initialS = pt.get_child("InitialState");
  auto InitialState = std::vector<int>(initialS.size());
  unsigned int counter(0);
  for (auto i : initialS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    unsigned int pos(counter++);
    boost::optional<unsigned int> opt_pos =
        i.second.get_optional<unsigned int>("<xmlattr>.PositionIndex");
    if (opt_pos)
      pos = opt_pos.get();
    InitialState.at(pos) = partP.getId();
  }

  auto finalS = pt.get_child("FinalState");
  auto FinalState = std::vector<int>(finalS.size());
  auto FinalStateEventPositionMapping =
      std::vector<unsigned int>(finalS.size());
  counter = 0;
  for (auto i : finalS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = partL->find(name)->second;
    unsigned int id = i.second.get<unsigned int>("<xmlattr>.Id");
    unsigned int pos(counter++);
    boost::optional<unsigned int> opt_pos =
        i.second.get_optional<unsigned int>("<xmlattr>.PositionIndex");
    if (opt_pos)
      pos = opt_pos.get();
    FinalState.at(pos) = partP.getId();
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

  return FourMomentum(px, py, pz, E);
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createIntensityFT(std::shared_ptr<PartList> partL,
                                       Kinematics &kin,
                                       const boost::property_tree::ptree &pt,
                                       std::string suffix) {
  LOG(TRACE) << "loading intensity...";

  std::string IntensityClass(pt.get<std::string>("<xmlattr>.Class"));

  if (IntensityClass == "IncoherentIntensity") {
    return createIncoherentIntensityFT(partL, kin, pt, suffix);
  } else if (IntensityClass == "CoherentIntensity") {
    return createCoherentIntensityFT(partL, kin, pt, suffix);
  } else if (IntensityClass == "StrengthIntensity") {
    return createStrengthIntensityFT(partL, kin, pt, suffix);
  } else if (IntensityClass == "NormalizedIntensity") {
    return createNormalizedIntensityFT(partL, kin, pt, suffix);
  } else {
    throw BadConfig(
        "IntensityBuilderXML::createIntensityFT() | Found unknown intensity " +
        IntensityClass);
  }
  return std::shared_ptr<ComPWA::FunctionTree::FunctionTree>(nullptr);
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createIncoherentIntensityFT(
    std::shared_ptr<PartList> partL, Kinematics &kin,
    const boost::property_tree::ptree &pt, std::string suffix) {

  LOG(TRACE) << "constructing IncoherentIntensity ...";
  auto name = pt.get<std::string>("<xmlattr>.Name");
  auto NodeName = "IncoherentIntensity(" + name + ")" + suffix;

  using namespace ComPWA::FunctionTree;

  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MDouble("", 0), std::make_shared<AddAll>(ParType::MDOUBLE));

  for (const auto &x : pt) {
    if (x.first == "Intensity") {
      tr->insertTree(createIntensityFT(partL, kin, x.second, suffix), NodeName);
    }
  }

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createCoherentIntensityFT(
    std::shared_ptr<PartList> partL, Kinematics &kin,
    const boost::property_tree::ptree &pt, std::string suffix) {

  LOG(TRACE) << "constructing CoherentIntensity ...";
  auto name = pt.get<std::string>("<xmlattr>.Name");

  auto NodeName = "CoherentIntensity(" + name + ")" + suffix;

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MDouble("", 0), std::make_shared<AbsSquare>(ParType::MDOUBLE));

  tr->createNode("SumOfAmplitudes" + suffix, MComplex("", 0),
                 std::make_shared<AddAll>(ParType::MCOMPLEX), NodeName);

  for (const auto &x : pt) {
    if (x.first == "Amplitude") {
      auto amp_tree = createAmplitudeFT(partL, kin, x.second, suffix);
      if (!amp_tree->sanityCheck())
        throw std::runtime_error(
            "CoherentIntensity::createFunctionTree(): tree "
            "didn't pass sanity check!");

      tr->insertTree(amp_tree, "SumOfAmplitudes" + suffix);
    }
  }

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createStrengthIntensityFT(
    std::shared_ptr<PartList> partL, Kinematics &kin,
    const boost::property_tree::ptree &pt, std::string suffix) {
  LOG(TRACE) << "creating StrengthIntensity ...";

  auto name = pt.get<std::string>("<xmlattr>.Name");

  std::shared_ptr<FitParameter> Strength(nullptr);
  boost::property_tree::ptree UndecoratedIntensityPT;
  for (const auto &x : pt) {
    if (x.first == "<xmlattr>")
      continue;
    if (x.first == "Parameter" &&
        x.second.get<std::string>("<xmlattr>.Type") == "Strength") {
      Strength = std::make_shared<FitParameter>(x.second);
    } else if (x.first == "Intensity") {
      UndecoratedIntensityPT = x.second;
    } else {
      LOG(WARNING) << "IntensityBuilderXML::createStrengthIntensity(): found "
                      "unknown tag "
                   << x.first;
    }
  }

  Strength = CurrentIntensityState.Parameters.addUniqueParameter(Strength);

  auto NodeName = "StrengthIntensity(" + name + ")" + suffix;

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MDouble("", 0), std::make_shared<MultAll>(ParType::MDOUBLE));

  tr->createLeaf("Strength", Strength, NodeName);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> x =
      createIntensityFT(partL, kin, UndecoratedIntensityPT, suffix);

  x->parameter();
  tr->insertTree(x, NodeName);

  if (!tr->sanityCheck())
    throw std::runtime_error(
        "StrengthIntensityDecorator::createFunctionTree() | "
        "Tree didn't pass sanity check!");

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createNormalizedIntensityFT(
    std::shared_ptr<PartList> partL, Kinematics &kin,
    const boost::property_tree::ptree &pt, std::string suffix) {
  LOG(TRACE) << "creating NormalizedIntensity ...";

  auto name = pt.get<std::string>("<xmlattr>.Name");

  boost::property_tree::ptree UnnormalizedPT;
  std::string IntegratorClassName("MCIntegrationStrategy");

  for (const auto &x : pt) {
    if (x.first == "<xmlattr>")
      continue;
    if (x.first == "Intensity") {
      UnnormalizedPT = x.second;
    } else if (x.first == "IntegrationStrategy") {
      auto OptionalIntegratorName =
          pt.get_optional<std::string>("<xmlattr>.IntegrationStrategy");
      if (OptionalIntegratorName.is_initialized()) {
        IntegratorClassName = OptionalIntegratorName.get();
      } else {
        LOG(INFO)
            << "IntensityBuilderXML::createNormalizedIntensityFT(): creating "
               "default IntegrationStrategy *MCIntegrationStrategy*";
      }
    } else {
      LOG(WARNING)
          << "IntensityBuilderXML::createNormalizedIntensityFT(): found "
             "unknown tag "
          << x.first;
    }
  }

  return normalizeIntensityFT(partL, kin, UnnormalizedPT, IntegratorClassName);
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::normalizeIntensityFT(
    std::shared_ptr<PartList> partL, Kinematics &kin,
    const boost::property_tree::ptree &UnnormalizedPT,
    std::string IntegratorClassName) {
  LOG(TRACE) << "creating Normalized FunctionTree ...";

  auto name = UnnormalizedPT.get<std::string>("<xmlattr>.Name");

  using namespace ComPWA::FunctionTree;

  auto NodeName = "NormalizedIntensity(" + name + ")";

  auto NormalizedFT = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MDouble("", 0), std::make_shared<MultAll>(ParType::MDOUBLE));

  // it is assumed that the ActiveData is already Data
  auto FTData = createIntensityFT(partL, kin, UnnormalizedPT, "");
  NormalizedFT->insertTree(FTData, NodeName);

  CurrentIntensityState.Data = CurrentIntensityState.ActiveData;
  CurrentIntensityState.ActiveData = CurrentIntensityState.PhspData;
  CurrentIntensityState.IsDataActive = false;
  auto FTPhspData = createIntensityFT(partL, kin, UnnormalizedPT, "_phsp");
  auto normtree =
      createIntegrationStrategyFT(FTPhspData, kin, IntegratorClassName);
  CurrentIntensityState.IsDataActive = true;
  CurrentIntensityState.PhspData = CurrentIntensityState.ActiveData;
  CurrentIntensityState.ActiveData = CurrentIntensityState.Data;

  NormalizedFT->insertTree(normtree, NodeName);

  return NormalizedFT;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createIntegrationStrategyFT(
    std::shared_ptr<ComPWA::FunctionTree::FunctionTree> UnnormalizedIntensity,
    Kinematics &kin, std::string IntegratorClassName) {
  LOG(TRACE) << "creating IntegrationStrategy ...";

  if (PhspSample.size() == 0)
    LOG(FATAL) << "IntensityBuilderXML::createIntegrationStrategyFT(): phsp "
                  "sample is not set!";

  using namespace ComPWA::FunctionTree;
  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> tr;

  if (IntegratorClassName == "MCIntegrationStrategy") {
    if (PhspSample.size() == 0)
      LOG(FATAL) << "IntensityBuilderXML::createIntegrationStrategyFT(): phsp "
                    "sample is not set!";

    // update the PhspData
    updateDataContainerState(kin);
    auto PhspDataSet = ComPWA::Data::convertEventsToDataSet(PhspSample, kin);
    ComPWA::FunctionTree::updateDataContainers(CurrentIntensityState.ActiveData,
                                               PhspDataSet.Data);

    if (!PhspWeights) {
      PhspWeights = MDouble("Weight", PhspDataSet.Weights);
    }
    double PhspWeightSum(std::accumulate(PhspWeights->values().begin(),
                                         PhspWeights->values().end(), 0.0));

    std::string NodeName =
        "Normalization(" + UnnormalizedIntensity->head()->name() + ")";
    tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
        NodeName, ValueFactory(ParType::DOUBLE),
        std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)));
    tr->createNode("Integral",
                   std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
                   NodeName);
    // normTree->createLeaf("PhspVolume", PhspVolume, "Integral");
    tr->createLeaf("InverseSampleWeights", 1.0 / PhspWeightSum, "Integral");
    tr->createNode("Sum",
                   std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                   "Integral");
    tr->createNode("WeightedIntensities",
                   std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)),
                   "Sum");

    if (PhspWeights)
      tr->createLeaf("EventWeight", PhspWeights, "WeightedIntensities");
    tr->insertTree(UnnormalizedIntensity, "WeightedIntensities");

    tr->parameter();
  } else {
    LOG(WARNING) << "IntensityBuilderXML::createIntegrationStrategyFT(): "
                    "IntegrationStrategy type "
                 << IntegratorClassName << " unknown!";
  }

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createAmplitudeFT(std::shared_ptr<PartList> partL,
                                       Kinematics &kin,
                                       const boost::property_tree::ptree &pt,
                                       std::string suffix) {
  auto ampclass = pt.get<std::string>("<xmlattr>.Class");

  if (ampclass == "HelicityDecay") {
    return createHelicityDecayFT(
        partL, *dynamic_cast<HelicityFormalism::HelicityKinematics *>(&kin), pt,
        suffix);
  } else if (ampclass == "CoefficientAmplitude") {
    return createCoefficientAmplitudeFT(partL, kin, pt, suffix);
  } else if (ampclass == "SequentialAmplitude") {
    return createSequentialAmplitudeFT(partL, kin, pt, suffix);
  } else if (ampclass == "NormalizedAmplitude") {
    return createNormalizedAmplitudeFT(partL, kin, pt);
  } else {
    throw BadConfig(
        "IntensityBuilderXML::createAmplitude(): Unknown amplitude " +
        ampclass);
  }
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createNormalizedAmplitudeFT(
    std::shared_ptr<PartList> partL, Kinematics &kin,
    const boost::property_tree::ptree &pt) {
  LOG(TRACE) << "creating NormalizedAmplitude ...";

  boost::property_tree::ptree UnnormalizedPT;
  std::string IntegratorClassName("MCIntegrationStrategy");

  for (const auto &x : pt) {
    if (x.first == "<xmlattr>")
      continue;
    if (x.first == "Amplitude") {
      UnnormalizedPT = x.second;
    } else if (x.first == "IntegrationStrategy") {
      auto OptionalIntegratorName =
          pt.get_optional<std::string>("<xmlattr>.IntegrationStrategy");
      if (OptionalIntegratorName.is_initialized()) {
        IntegratorClassName = OptionalIntegratorName.get();
      } else {
        LOG(INFO)
            << "IntensityBuilderXML::createNormalizedAmplitudeFT(): creating "
               "default IntegrationStrategy *MCIntegrationStrategy*";
      }
    } else {
      LOG(WARNING)
          << "IntensityBuilderXML::createNormalizedAmplitudeFT(): found "
             "unknown tag "
          << x.first;
    }
  }

  auto FTData = createAmplitudeFT(partL, kin, UnnormalizedPT, "");

  auto name = FTData->head()->name();

  auto NodeName = "Normalized(" + name + ")";

  using namespace ComPWA::FunctionTree;

  auto NormalizedFT = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MComplex("", 0), std::make_shared<MultAll>(ParType::MCOMPLEX));

  NormalizedFT->insertTree(FTData, NodeName);

  // if we are using the Data as ActiveData, then swap
  if (CurrentIntensityState.IsDataActive) {
    CurrentIntensityState.Data = CurrentIntensityState.ActiveData;
    CurrentIntensityState.ActiveData = CurrentIntensityState.PhspData;
  }

  auto FTPhspData = createAmplitudeFT(partL, kin, UnnormalizedPT, "_phsp");
  // this phspdata function tree has to be made into a double valued function
  auto FTPhspDataAbsSquared =
      std::make_shared<ComPWA::FunctionTree::FunctionTree>(
          FTPhspData->head()->name() + "_AbsSquared", MDouble("", 0),
          std::make_shared<AbsSquare>(ParType::MDOUBLE));
  FTPhspDataAbsSquared->insertTree(FTPhspData,
                                   FTPhspData->head()->name() + "_AbsSquared");

  auto normtreesquared = createIntegrationStrategyFT(FTPhspDataAbsSquared, kin,
                                                     IntegratorClassName);

  auto normtree = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      normtreesquared->head()->name() + "_sqrt", ValueFactory(ParType::DOUBLE),
      std::make_shared<SquareRoot>(ParType::DOUBLE));

  normtree->insertTree(normtreesquared,
                       normtreesquared->head()->name() + "_sqrt");

  if (CurrentIntensityState.IsDataActive) {
    CurrentIntensityState.PhspData = CurrentIntensityState.ActiveData;
    CurrentIntensityState.ActiveData = CurrentIntensityState.Data;
  }
  NormalizedFT->insertTree(normtree, NodeName);

  return NormalizedFT;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createCoefficientAmplitudeFT(
    std::shared_ptr<PartList> partL, Kinematics &kin,
    const boost::property_tree::ptree &pt, std::string suffix) {
  LOG(TRACE) << "constructing CoefficientAmplitude ...";

  std::shared_ptr<FitParameter> Magnitude(nullptr);
  std::shared_ptr<FitParameter> Phase(nullptr);
  boost::property_tree::ptree AmpPT;
  for (const auto &v : pt) {
    if (v.first == "Parameter") {
      if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude")
        Magnitude = std::make_shared<FitParameter>(v.second);
      if (v.second.get<std::string>("<xmlattr>.Type") == "Phase")
        Phase = std::make_shared<FitParameter>(v.second);
    } else if (v.first == "Amplitude") {
      AmpPT = v.second;
    }
  }

  if (!Magnitude)
    throw BadParameter(
        "IntensityBuilderXML::createCoefficientAmplitude() | No magnitude "
        "parameter found.");
  if (!Phase)
    throw BadParameter("IntensityBuilderXML::createCoefficientAmplitude() | No "
                       "phase parameter found.");

  auto amp_ft = createAmplitudeFT(partL, kin, AmpPT, suffix);

  std::string ampname = amp_ft->head()->name();

  Magnitude = CurrentIntensityState.Parameters.addUniqueParameter(Magnitude);
  Phase = CurrentIntensityState.Parameters.addUniqueParameter(Phase);

  std::string nodeName = "CoefficientAmplitude(" + ampname + ")" + suffix;

  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      nodeName, ComPWA::FunctionTree::MComplex("", 0),
      std::make_shared<ComPWA::FunctionTree::MultAll>(
          ComPWA::FunctionTree::ParType::MCOMPLEX));
  tr->createNode(
      "Strength",
      std::make_shared<ComPWA::FunctionTree::Value<std::complex<double>>>(),
      std::make_shared<ComPWA::FunctionTree::Complexify>(
          ComPWA::FunctionTree::ParType::COMPLEX),
      nodeName);
  tr->createLeaf("Magnitude", Magnitude, "Strength");
  tr->createLeaf("Phase", Phase, "Strength");

  tr->insertTree(amp_ft, nodeName);
  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createSequentialAmplitudeFT(
    std::shared_ptr<PartList> partL, Kinematics &kin,
    const boost::property_tree::ptree &pt, std::string suffix) {
  LOG(TRACE) << "constructing SequentialAmplitude ...";

  std::string ampname;

  std::vector<std::shared_ptr<ComPWA::FunctionTree::FunctionTree>> Amplitudes;
  auto PreFactor = std::complex<double>(1, 0);
  for (const auto &v : pt) {
    if (v.first == "Amplitude") {
      std::shared_ptr<ComPWA::FunctionTree::FunctionTree> AmpTree =
          createAmplitudeFT(partL, kin, v.second, suffix);
      if (!AmpTree->sanityCheck())
        throw std::runtime_error(
            "SequentialAmplitude::createFunctionTree : tree "
            "didn't pass sanity check!");
      AmpTree->parameter();
      ampname += AmpTree->head()->name();
      Amplitudes.push_back(AmpTree);
    } else if (v.first == "PreFactor") {
      boost::optional<double> optr =
          v.second.get_optional<double>("<xmlattr>.Magnitude");
      if (optr.is_initialized()) {
        double r(optr.value());
        if (r < 0.0)
          throw BadConfig("IntensityBuilderXML::createSequentialAmplitude(): "
                          "PreFactor Magnitude below zero!");
        double p(0.0);
        boost::optional<double> optp =
            v.second.get_optional<double>("<xmlattr>.Phase");
        if (optp.is_initialized())
          p = optp.value();
        PreFactor = std::polar(r, p);
      } else {
        double real = v.second.get<double>("<xmlattr>.Real");
        double im(0.0);
        boost::optional<double> optim =
            v.second.get_optional<double>("<xmlattr>.Imaginary");
        if (optim.is_initialized())
          im = optim.value();
        PreFactor = std::complex<double>(real, im);
      }
    } else if (v.first != "<xmlattr>") {
      throw BadConfig("SequentialAmplitude::createSequentialAmplitude() | "
                      "Unknown tag " +
                      v.first + "!");
    }
  }

  auto NodeName = "SequentialPartialAmplitude(" + ampname + ")" + suffix;

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MComplex("", 0), std::make_shared<MultAll>(ParType::MCOMPLEX));

  for (auto x : Amplitudes) {
    tr->insertTree(x, NodeName);
  }
  tr->createLeaf("Prefactor", PreFactor, NodeName);

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createHelicityDecayFT(
    std::shared_ptr<PartList> partL, HelicityFormalism::HelicityKinematics &kin,
    const boost::property_tree::ptree &pt, std::string suffix) {
  LOG(TRACE) << "IntensityBuilderXML::createHelicityDecayFT(): ";
  unsigned int SubSystemIndex(kin.addSubSystem(SubSystem(pt)));
  unsigned int DataPosition = 3 * SubSystemIndex;

  updateDataContainerState(kin);

  auto ampname = pt.get<std::string>("<xmlattr>.Name") + "_" +
                 std::to_string(DataPosition) + ";";

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  auto partItr = partL->find(name);
  if (partItr == partL->end())
    throw std::runtime_error(
        "IntensityBuilderXML::createHelicityDecayFT(): Particle " + name +
        " not found in list!");
  ComPWA::Spin J = partItr->second.getSpinQuantumNumber("Spin");
  ComPWA::Spin mu(pt.get<double>("DecayParticle.<xmlattr>.Helicity"));
  // if the node OrbitalAngularMomentum does not exist, set it to spin J as
  // default value
  ComPWA::Spin orbitL(pt.get<double>(
      "DecayParticle.<xmlattr>.OrbitalAngularMomentum", (double)J));

  auto partProp = partItr->second;
  auto Mass = std::make_shared<FunctionTree::FitParameter>(
      partProp.getMass().Name, partProp.getMass().Value,
      partProp.getMass().Error.first);
  Mass->fixParameter(partProp.getMass().IsFixed);
  Mass = CurrentIntensityState.Parameters.addUniqueParameter(Mass);

  double PreFactor(1.0);
  const auto &canoSum = pt.get_child_optional("CanonicalSum");
  if (canoSum) {
    const auto &sumTree = canoSum.get();
    orbitL = sumTree.get<double>("<xmlattr>.L");
    double coef = std::sqrt((2 * (double)orbitL + 1) / (2 * (double)J + 1));
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
        "IntensityBuilderXML::createHelicityDecayFT(): Expect exactly two "
        "decay products (" +
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

  auto parMass1 = std::make_shared<FitParameter>(
      partL->find(DecayProducts.first)->second.getMass().Name,
      partL->find(DecayProducts.first)->second.getMass().Value);
  auto parMass2 = std::make_shared<FitParameter>(
      partL->find(DecayProducts.second)->second.getMass().Name,
      partL->find(DecayProducts.second)->second.getMass().Value);
  parMass1->fixParameter(
      partL->find(DecayProducts.first)->second.getMass().IsFixed);
  parMass2->fixParameter(
      partL->find(DecayProducts.second)->second.getMass().IsFixed);
  parMass1 = CurrentIntensityState.Parameters.addUniqueParameter(parMass1);
  parMass2 = CurrentIntensityState.Parameters.addUniqueParameter(parMass2);
  auto decayInfo = partProp.getDecayInfo();
  Dynamics::FormFactorType ffType = Dynamics::FormFactorType::noFormFactor;
  std::shared_ptr<FitParameter> parRadius;
  std::shared_ptr<FitParameter> Width;
  for (const auto &node : decayInfo.get_child("")) {
    if (node.first == "FormFactor") {
      auto FFTypeInt = node.second.get<int>("<xmlattr>.Type");
      if (Dynamics::FormFactorType::BlattWeisskopf == FFTypeInt)
        ffType = Dynamics::FormFactorType::BlattWeisskopf;
      else if (Dynamics::FormFactorType::CrystalBarrel == FFTypeInt)
        ffType = Dynamics::FormFactorType::CrystalBarrel;
    } else if (node.first == "Parameter") {
      std::string parType = node.second.get<std::string>("<xmlattr>.Type");
      if (parType == "Width") {
        Width = std::make_shared<FitParameter>(node.second);
        Width = CurrentIntensityState.Parameters.addUniqueParameter(Width);
      } else if (parType == "MesonRadius") {
        parRadius = std::make_shared<FitParameter>(node.second);
        parRadius =
            CurrentIntensityState.Parameters.addUniqueParameter(parRadius);
      }
    }
  }

  std::string decayType = partProp.getDecayType();

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> DynamicFunctionFT(
      nullptr);

  if (decayType == "stable") {
    throw std::runtime_error(
        "IntensityBuilderXML::createHelicityDecayFT(): Stable particle is "
        "given as mother particle of a decay. Makes no sense!");
  } else if (decayType == "relativisticBreitWigner") {
    using namespace ComPWA::Physics::Dynamics::RelativisticBreitWigner;
    InputInfo RBW;
    RBW.Mass = Mass;
    RBW.Width = Width;
    RBW.MesonRadius = parRadius;
    RBW.DaughterMasses = std::make_pair(parMass1, parMass2);
    RBW.FFType = ffType;
    RBW.L = (unsigned int)orbitL;
    DynamicFunctionFT = createFunctionTree(
        RBW, CurrentIntensityState.ActiveData, DataPosition, suffix);
  } else if (decayType == "flatte") {
    using namespace ComPWA::Physics::Dynamics::Flatte;
    InputInfo FlatteInfo;
    FlatteInfo.Mass = Mass;
    FlatteInfo.Width = Width;
    FlatteInfo.MesonRadius = parRadius;
    FlatteInfo.DaughterMasses = std::make_pair(parMass1, parMass2);
    FlatteInfo.FFType = ffType;
    FlatteInfo.L = (unsigned int)orbitL;
    // Read parameters
    for (const auto &v : decayInfo.get_child("")) {
      if (v.first != "Parameter")
        continue;
      std::string type = v.second.get<std::string>("<xmlattr>.Type");
      if (type == "Coupling") {
        FlatteInfo.Couplings.push_back(Dynamics::Coupling(partL, v.second));
      }
    }
    DynamicFunctionFT = createFunctionTree(
        FlatteInfo, CurrentIntensityState.ActiveData, DataPosition, suffix);
  } else if (decayType == "voigt") {
    using namespace ComPWA::Physics::Dynamics::Voigtian;
    InputInfo VoigtInfo;
    VoigtInfo.Mass = Mass;
    VoigtInfo.Width = Width;
    VoigtInfo.MesonRadius = parRadius;
    VoigtInfo.DaughterMasses = std::make_pair(parMass1, parMass2);
    VoigtInfo.FFType = ffType;
    VoigtInfo.L = (unsigned int)orbitL;
    VoigtInfo.Sigma = decayInfo.get<double>("Resolution.<xmlattr>.Sigma");
    DynamicFunctionFT = createFunctionTree(
        VoigtInfo, CurrentIntensityState.ActiveData, DataPosition, suffix);
  } else if (decayType == "virtual" || decayType == "nonResonant") {
    DynamicFunctionFT = Dynamics::NonResonant::createFunctionTree(
        CurrentIntensityState.ActiveData, DataPosition, suffix);
  } else {
    throw std::runtime_error("HelicityDecay::Factory() | Unknown decay type " +
                             decayType + "!");
  }

  auto AngularFunction =
      ComPWA::Physics::HelicityFormalism::WignerD::createFunctionTree(
          J, mu, DecayHelicities.first - DecayHelicities.second,
          CurrentIntensityState.ActiveData, DataPosition + 1, DataPosition + 2);

  std::string nodeName = "PartialAmplitude(" + ampname + ")" + suffix;

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      nodeName, MComplex("", 0), std::make_shared<MultAll>(ParType::MCOMPLEX));
  tr->createLeaf("PreFactor", PreFactor, nodeName);
  tr->insertTree(AngularFunction, nodeName);
  tr->insertTree(DynamicFunctionFT, nodeName);

  // set production formfactor
  if (ffType != ComPWA::Physics::Dynamics::FormFactorType::noFormFactor &&
      ((unsigned int)orbitL > 0)) {
    if (parRadius == nullptr) {
      throw std::runtime_error("IntensityBuilderXML::createHelicityDecayFT(): "
                               "No MesonRadius given for amplitude " +
                               name +
                               "! It is needed to calculate the form factor!");
    }

    std::shared_ptr<ComPWA::FunctionTree::FunctionTree> ProductionFormFactorFT =
        Dynamics::createFunctionTree(
            name, parMass1, parMass2, parRadius, orbitL, ffType,
            CurrentIntensityState.ActiveData, DataPosition, suffix);

    tr->insertTree(ProductionFormFactorFT, nodeName);
  }

  tr->parameter();

  return tr;
}

void IntensityBuilderXML::updateDataContainerState(const Kinematics &kin) {
  LOG(TRACE) << "updating data container...";
  LOG(TRACE) << "currently there are "
             << CurrentIntensityState.ActiveData.mDoubleValues().size()
             << " entries."
             << " kinematics has " << kin.getKinematicVariableNames().size()
             << " entries.";

  for (size_t i = CurrentIntensityState.ActiveData.mDoubleValues().size();
       i < kin.getKinematicVariableNames().size(); ++i) {
    std::vector<double> temp;
    CurrentIntensityState.ActiveData.addValue(
        std::make_shared<ComPWA::FunctionTree::Value<std::vector<double>>>(
            kin.getKinematicVariableNames()[i], temp));
  }
}

} // namespace Physics
} // namespace ComPWA
