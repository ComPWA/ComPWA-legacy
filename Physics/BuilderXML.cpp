// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "BuilderXML.hpp"

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

IntensityBuilderXML::IntensityBuilderXML(
    ParticleList PartList_, Kinematics &Kin,
    const boost::property_tree::ptree &ModelTree_,
    const std::vector<Event> &TruePhspSample_,
    const std::vector<Event> &RecoPhspSample_)
    : PartList(PartList_), Kinematic(Kin), ModelTree(ModelTree_),
      TruePhspSample(TruePhspSample_),
      RecoPhspSample(
          [&TruePhspSample_,
           &RecoPhspSample_]() -> const std::vector<ComPWA::Event> & {
            if (TruePhspSample_.size() > 0 && RecoPhspSample_.size() == 0)
              return TruePhspSample_;
            else {
              return RecoPhspSample_;
            }
          }()) {}

ComPWA::FunctionTree::FunctionTreeIntensity
IntensityBuilderXML::createIntensity() {
  LOG(TRACE) << "loading intensity...";

  // BlankState
  Parameters = FunctionTree::ParameterList();
  PhspData = DataContainer();
  PhspRecoData = DataContainer();

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> FT =
      createIntensityFT(ModelTree, Data.Data, "");

  return {FT, Parameters, Data.Data};
}

std::vector<ComPWA::Tools::IntensityComponent>
IntensityBuilderXML::createIntensityComponents(
    std::vector<std::vector<std::string>> ComponentList) {
  LOG(TRACE) << "Creating intensity components...";

  if (UniqueComponentFTMapping.size() == 0) {
    LOG(INFO) << "Components map is empty. Creating full Intensity first.";
    createIntensity();
  }

  // then put these FT together according to the string vectors
  std::vector<ComPWA::Tools::IntensityComponent> IntensityComponents;

  for (auto const &Component : ComponentList) {
    std::string ComponentName;
    std::map<std::string,
             std::pair<std::string,
                       std::shared_ptr<ComPWA::FunctionTree::FunctionTree>>>
        NewUniqueComponentFTMapping;
    std::string Type("");
    std::vector<std::shared_ptr<ComPWA::FunctionTree::FunctionTree>> FTList;
    for (auto const &x : Component) {
      auto FindResult = UniqueComponentFTMapping.find(x);
      if (UniqueComponentFTMapping.end() != FindResult) {
        ComponentName += "_" + x;
        if (Type == "") {
          Type = FindResult->second.first;
        } else {
          if (Type != FindResult->second.first) {
            LOG(ERROR) << "Component " << x
                       << " incompatible with previous type " << Type
                       << " Skipping...";
            continue;
          }
        }
        NewUniqueComponentFTMapping.insert(*FindResult);
        FTList.push_back(FindResult->second.second);
      } else {
        LOG(ERROR) << "Component " << x << " not found! Skipping...";
      }
    }
    if (FTList.size() == 0)
      continue;

    ComponentName.erase(0, 1);

    LOG(INFO) << "Building component " << ComponentName;

    // empty component were already removed so [0] is safe
    if (Type == "Amplitude") {
      IntensityComponents.push_back(std::make_pair(
          ComponentName,
          std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensity>(
              createCoherentIntensityFT(ComponentName, FTList, ""), Parameters,
              Data.Data)));
      LOG(INFO) << "as a CoherentIntensity";
    } else {
      IntensityComponents.push_back(std::make_pair(
          ComponentName,
          std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensity>(
              createIncoherentIntensityFT(ComponentName, FTList, ""),
              Parameters, Data.Data)));
      LOG(INFO) << "as a IncoherentIntensity";
    }
  }
  return IntensityComponents;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createIntensityFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample, std::string suffix) {
  LOG(TRACE) << "loading intensity...";

  std::string IntensityClass(pt.get<std::string>("<xmlattr>.Class"));

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> FT(nullptr);

  if (IntensityClass == "IncoherentIntensity") {
    FT = createIncoherentIntensityFT(pt, DataSample, suffix);
  } else if (IntensityClass == "CoherentIntensity") {
    FT = createCoherentIntensityFT(pt, DataSample, suffix);
  } else if (IntensityClass == "StrengthIntensity") {
    FT = createStrengthIntensityFT(pt, DataSample, suffix);
  } else if (IntensityClass == "NormalizedIntensity") {
    FT = createNormalizedIntensityFT(pt, DataSample, suffix);
  } else {
    throw BadConfig("IntensityBuilderXML::createIntensityFT() | Found "
                    "unknown intensity " +
                    IntensityClass);
  }

  if ("" == suffix) {
    auto Name = pt.get<std::string>("<xmlattr>.Name");
    addFunctionTreeComponent(Name, "Intensity", FT);
  }

  return FT;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createIncoherentIntensityFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample, std::string suffix) {
  LOG(TRACE) << "constructing IncoherentIntensity ...";
  auto name = pt.get<std::string>("<xmlattr>.Name");

  std::vector<std::shared_ptr<ComPWA::FunctionTree::FunctionTree>> intens;
  for (const auto &x : pt) {
    if (x.first == "Intensity") {
      intens.push_back(createIntensityFT(x.second, DataSample, suffix));
    }
  }

  return createIncoherentIntensityFT(name, intens, suffix);
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createIncoherentIntensityFT(
    std::string Name,
    std::vector<std::shared_ptr<ComPWA::FunctionTree::FunctionTree>>
        Intensities,
    std::string suffix) {
  LOG(TRACE) << "constructing IncoherentIntensity ...";
  auto NodeName = "IncoherentIntensity(" + Name + ")" + suffix;

  using namespace ComPWA::FunctionTree;

  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MDouble("", 0), std::make_shared<AddAll>(ParType::MDOUBLE));

  for (auto x : Intensities) {
    tr->insertTree(x, NodeName);
  }

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createCoherentIntensityFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample, std::string suffix) {
  LOG(TRACE) << "constructing CoherentIntensity ...";
  auto name = pt.get<std::string>("<xmlattr>.Name");

  std::vector<std::shared_ptr<ComPWA::FunctionTree::FunctionTree>> amps;
  for (const auto &x : pt) {
    if (x.first == "Amplitude") {
      amps.push_back(createAmplitudeFT(x.second, DataSample, suffix));
    }
  }

  return createCoherentIntensityFT(name, amps, suffix);
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createCoherentIntensityFT(
    std::string Name,
    std::vector<std::shared_ptr<ComPWA::FunctionTree::FunctionTree>> Amplitudes,
    std::string suffix) {
  LOG(TRACE) << "constructing CoherentIntensity ...";

  auto NodeName = "CoherentIntensity(" + Name + ")" + suffix;

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MDouble("", 0), std::make_shared<AbsSquare>(ParType::MDOUBLE));

  tr->createNode("SumOfAmplitudes" + suffix, MComplex("", 0),
                 std::make_shared<AddAll>(ParType::MCOMPLEX), NodeName);

  for (auto x : Amplitudes) {
    tr->insertTree(x, "SumOfAmplitudes" + suffix);
  }

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createStrengthIntensityFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample, std::string suffix) {
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

  Strength = Parameters.addUniqueParameter(Strength);

  auto NodeName = "StrengthIntensity(" + name + ")" + suffix;

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MDouble("", 0), std::make_shared<MultAll>(ParType::MDOUBLE));

  tr->createLeaf("Strength", Strength, NodeName);

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> x =
      createIntensityFT(UndecoratedIntensityPT, DataSample, suffix);

  x->parameter();
  tr->insertTree(x, NodeName);

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createNormalizedIntensityFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample, std::string suffix) {
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

  return normalizeIntensityFT(UnnormalizedPT, DataSample, IntegratorClassName,
                              suffix);
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::normalizeIntensityFT(
    const boost::property_tree::ptree &UnnormalizedPT,
    const ComPWA::FunctionTree::ParameterList &DataSample,
    std::string IntegratorClassName, std::string suffix) {
  LOG(TRACE) << "creating Normalized FunctionTree ...";

  if (RecoPhspSample.size() == 0)
    LOG(FATAL) << "IntensityBuilderXML::normalizeIntensityFT(): "
                  "reco phsp sample is not set!";
  updateDataContainerWeights(PhspRecoData, RecoPhspSample);

  auto name = UnnormalizedPT.get<std::string>("<xmlattr>.Name");

  using namespace ComPWA::FunctionTree;

  auto NodeName = "NormalizedIntensity(" + name + ")";

  auto NormalizedFT = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MDouble("", 0), std::make_shared<MultAll>(ParType::MDOUBLE));

  auto FTData = createIntensityFT(UnnormalizedPT, DataSample, suffix);
  NormalizedFT->insertTree(FTData, NodeName);

  auto FTPhspData =
      createIntensityFT(UnnormalizedPT, PhspRecoData.Data, "_phsp_rec");
  auto normtree =
      createIntegrationStrategyFT(FTPhspData, PhspRecoData.Weights,
                                  PhspRecoData.WeightSum, IntegratorClassName);

  NormalizedFT->insertTree(normtree, NodeName);

  return NormalizedFT;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createIntegrationStrategyFT(
    std::shared_ptr<ComPWA::FunctionTree::FunctionTree> UnnormalizedIntensity,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        PhspWeights,
    double PhspWeightSum, std::string IntegratorClassName) {
  LOG(TRACE) << "creating IntegrationStrategy ...";

  using namespace ComPWA::FunctionTree;
  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> tr;

  if (IntegratorClassName == "MCIntegrationStrategy") {
    // update the PhspData container
    updateDataContainerState();
    updateDataContainerContent();

    std::string NodeName =
        "Normalization(" + UnnormalizedIntensity->Head->name() + ")";
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
IntensityBuilderXML::createAmplitudeFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample, std::string suffix) {
  auto ampclass = pt.get<std::string>("<xmlattr>.Class");

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree> FT(nullptr);

  if (ampclass == "HelicityDecay") {
    FT = createHelicityDecayFT(pt, DataSample, suffix);
  } else if (ampclass == "CoefficientAmplitude") {
    FT = createCoefficientAmplitudeFT(pt, DataSample, suffix);
  } else if (ampclass == "SequentialAmplitude") {
    FT = createSequentialAmplitudeFT(pt, DataSample, suffix);
  } else if (ampclass == "NormalizedAmplitude") {
    FT = createNormalizedAmplitudeFT(pt, DataSample, suffix);
  } else {
    throw BadConfig(
        "IntensityBuilderXML::createAmplitude(): Unknown amplitude " +
        ampclass);
  }

  if ("" == suffix) {
    auto Name = pt.get<std::string>("<xmlattr>.Name");
    addFunctionTreeComponent(Name, "Amplitude", FT);
  }

  return FT;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createNormalizedAmplitudeFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample, std::string suffix) {
  LOG(TRACE) << "creating NormalizedAmplitude ...";

  if (TruePhspSample.size() == 0)
    LOG(FATAL) << "IntensityBuilderXML::createNormalizedAmplitudeFT(): "
                  "true phsp sample is not set!";
  updateDataContainerWeights(PhspData, TruePhspSample);

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

  auto FTData = createAmplitudeFT(UnnormalizedPT, DataSample, suffix);

  auto name = FTData->Head->name();

  auto NodeName = "Normalized(" + name + ")";

  using namespace ComPWA::FunctionTree;

  auto NormalizedFT = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MComplex("", 0), std::make_shared<MultAll>(ParType::MCOMPLEX));

  NormalizedFT->insertTree(FTData, NodeName);

  auto FTPhspData = createAmplitudeFT(UnnormalizedPT, PhspData.Data, "_phsp");
  // this phspdata function tree has to be made into a double valued function
  auto FTPhspDataAbsSquared =
      std::make_shared<ComPWA::FunctionTree::FunctionTree>(
          FTPhspData->Head->name() + "_AbsSquared", MDouble("", 0),
          std::make_shared<AbsSquare>(ParType::MDOUBLE));
  FTPhspDataAbsSquared->insertTree(FTPhspData,
                                   FTPhspData->Head->name() + "_AbsSquared");

  auto normtreesquared =
      createIntegrationStrategyFT(FTPhspDataAbsSquared, PhspData.Weights,
                                  PhspData.WeightSum, IntegratorClassName);

  auto normtree = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      normtreesquared->Head->name() + "_sqrt", ValueFactory(ParType::DOUBLE),
      std::make_shared<SquareRoot>(ParType::DOUBLE));

  normtree->insertTree(normtreesquared,
                       normtreesquared->Head->name() + "_sqrt");

  NormalizedFT->insertTree(normtree, NodeName);

  return NormalizedFT;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
IntensityBuilderXML::createCoefficientAmplitudeFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample, std::string suffix) {
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

  auto amp_ft = createAmplitudeFT(AmpPT, DataSample, suffix);

  std::string ampname = amp_ft->Head->name();

  Magnitude = Parameters.addUniqueParameter(Magnitude);
  Phase = Parameters.addUniqueParameter(Phase);

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
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample, std::string suffix) {
  LOG(TRACE) << "constructing SequentialAmplitude ...";

  std::string ampname;

  std::vector<std::shared_ptr<ComPWA::FunctionTree::FunctionTree>> Amplitudes;
  auto PreFactor = std::complex<double>(1, 0);
  for (const auto &v : pt) {
    if (v.first == "Amplitude") {
      std::shared_ptr<ComPWA::FunctionTree::FunctionTree> AmpTree =
          createAmplitudeFT(v.second, DataSample, suffix);

      AmpTree->parameter();
      ampname += AmpTree->Head->name();
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
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample, std::string suffix) {
  LOG(TRACE) << "IntensityBuilderXML::createHelicityDecayFT(): ";
  auto &kin = dynamic_cast<HelicityFormalism::HelicityKinematics &>(Kinematic);
  unsigned int SubSystemIndex(kin.addSubSystem(SubSystem(pt)));
  unsigned int DataPosition = 3 * SubSystemIndex;

  updateDataContainerState();

  auto ampname = pt.get<std::string>("<xmlattr>.Name") + "_" +
                 std::to_string(DataPosition) + ";";

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  auto partProp = ComPWA::findParticle(PartList, name);

  double J = partProp.getQuantumNumber<double>("Spin");
  double mu(pt.get<double>("DecayParticle.<xmlattr>.Helicity"));
  // if the node OrbitalAngularMomentum does not exist, set it to spin J as
  // default value
  unsigned int orbitL(J);

  auto Mass = std::make_shared<FunctionTree::FitParameter>(
      partProp.getMass().Name, partProp.getMass().Value,
      partProp.getMass().Error.first);
  Mass->fixParameter(partProp.getMass().IsFixed);
  Mass = Parameters.addUniqueParameter(Mass);

  double PreFactor(1.0);
  const auto &canoSum = pt.get_child_optional("CanonicalSum");
  if (canoSum) {
    const auto &sumTree = canoSum.get();
    orbitL = sumTree.get<unsigned int>("<xmlattr>.L");
    double coef = std::sqrt((2.0 * orbitL + 1) / (2 * J + 1));
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
  std::pair<double, double> DecayHelicities;

  auto p = decayProducts.begin();
  DecayProducts.first = p->second.get<std::string>("<xmlattr>.Name");
  DecayHelicities.first = p->second.get<double>("<xmlattr>.Helicity");
  ++p;
  DecayProducts.second = p->second.get<std::string>("<xmlattr>.Name");
  DecayHelicities.second = p->second.get<double>("<xmlattr>.Helicity");

  auto parMass1 = std::make_shared<FitParameter>(
      ComPWA::findParticle(PartList, DecayProducts.first).getMass().Name,
      ComPWA::findParticle(PartList, DecayProducts.first).getMass().Value);
  auto parMass2 = std::make_shared<FitParameter>(
      ComPWA::findParticle(PartList, DecayProducts.second).getMass().Name,
      ComPWA::findParticle(PartList, DecayProducts.second).getMass().Value);
  parMass1->fixParameter(
      ComPWA::findParticle(PartList, DecayProducts.first).getMass().IsFixed);
  parMass2->fixParameter(
      ComPWA::findParticle(PartList, DecayProducts.second).getMass().IsFixed);
  parMass1 = Parameters.addUniqueParameter(parMass1);
  parMass2 = Parameters.addUniqueParameter(parMass2);
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
        Width = Parameters.addUniqueParameter(Width);
      } else if (parType == "MesonRadius") {
        parRadius = std::make_shared<FitParameter>(node.second);
        parRadius = Parameters.addUniqueParameter(parRadius);
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
    DynamicFunctionFT =
        createFunctionTree(RBW, DataSample, DataPosition, suffix);
  } else if (decayType == "flatte") {
    ComPWA::Physics::Dynamics::Flatte::InputInfo FlatteInfo;
    FlatteInfo.Mass = Mass;
    FlatteInfo.MesonRadius = parRadius;
    FlatteInfo.DaughterMasses = std::make_pair(parMass1, parMass2);
    FlatteInfo.FFType = ffType;
    FlatteInfo.L = (unsigned int)orbitL;
    std::vector<Dynamics::Coupling> couplings;
    // Read parameters
    for (const auto &v : decayInfo.get_child("")) {
      if (v.first != "Parameter")
        continue;
      std::string type = v.second.get<std::string>("<xmlattr>.Type");
      if (type == "Coupling") {
        auto c = Dynamics::Coupling(PartList, v.second);
        c.G = Parameters.addUniqueParameter(c.G);

        if ((c.MassA->value() == parMass1->value() &&
             c.MassB->value() == parMass2->value()) ||
            (c.MassB->value() == parMass1->value() &&
             c.MassA->value() == parMass2->value())) {
          FlatteInfo.G = c.G;
        } else {
          c.MassA = Parameters.addUniqueParameter(c.MassA);
          c.MassB = Parameters.addUniqueParameter(c.MassB);
          couplings.push_back(c);
        }
      }
    }
    FlatteInfo.HiddenCouplings = couplings;
    DynamicFunctionFT =
        createFunctionTree(FlatteInfo, DataSample, DataPosition, suffix);
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
    DynamicFunctionFT =
        createFunctionTree(VoigtInfo, DataSample, DataPosition, suffix);
  } else if (decayType == "virtual" || decayType == "nonResonant") {
    DynamicFunctionFT = Dynamics::NonResonant::createFunctionTree(
        DataSample, DataPosition, suffix);
  } else {
    throw std::runtime_error("HelicityDecay::Factory() | Unknown decay type " +
                             decayType + "!");
  }

  auto AngularFunction =
      ComPWA::Physics::HelicityFormalism::WignerD::createFunctionTree(
          J, mu, DecayHelicities.first - DecayHelicities.second, DataSample,
          DataPosition + 1, DataPosition + 2);

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
        Dynamics::createFunctionTree(name, parMass1, parMass2, parRadius,
                                     orbitL, ffType, DataSample, DataPosition,
                                     suffix);

    tr->insertTree(ProductionFormFactorFT, nodeName);
  }

  tr->parameter();

  return tr;
}

void updateDataContainerState(ComPWA::FunctionTree::ParameterList &DataSample,
                              const Kinematics &Kin) {
  LOG(TRACE) << "currently there are " << DataSample.mDoubleValues().size()
             << " entries."
             << " kinematics has " << Kin.getKinematicVariableNames().size()
             << " entries.";

  for (size_t i = DataSample.mDoubleValues().size();
       i < Kin.getKinematicVariableNames().size(); ++i) {
    std::vector<double> temp;
    DataSample.addValue(
        std::make_shared<ComPWA::FunctionTree::Value<std::vector<double>>>(
            Kin.getKinematicVariableNames()[i], temp));
  }
}

void IntensityBuilderXML::updateDataContainerState() {
  LOG(TRACE) << "updating data containers...";
  Physics::updateDataContainerState(Data.Data, Kinematic);
  Physics::updateDataContainerState(PhspData.Data, Kinematic);
  Physics::updateDataContainerState(PhspRecoData.Data, Kinematic);
}

void updateDataContainerContent(ComPWA::FunctionTree::ParameterList &DataList,
                                const std::vector<ComPWA::Event> &DataSample,
                                const Kinematics &Kin) {
  auto DataSet = ComPWA::Data::convertEventsToDataSet(DataSample, Kin);

  // just loop over the vectors and fill in the data
  if (DataList.mDoubleValues().size() > DataSet.Data.size()) {
    std::stringstream ss;
    ss << "IntensityBuilderXML::updateDataContainerContent(): given data "
          "container does not have enough variables! (required: "
       << DataList.mDoubleValues().size() << ", given: " << DataSet.Data.size()
       << ")";
    throw std::out_of_range(ss.str());
  }
  for (size_t i = 0; i < DataList.mDoubleValues().size(); ++i) {
    if (DataList.mDoubleValue(i)->values().size() == 0) {
      DataList.mDoubleValue(i)->setValue(DataSet.Data[i]);
    }
  }
}

void IntensityBuilderXML::updateDataContainerContent() {
  Physics::updateDataContainerContent(PhspData.Data, TruePhspSample, Kinematic);
  Physics::updateDataContainerContent(PhspRecoData.Data, RecoPhspSample,
                                      Kinematic);
}

void IntensityBuilderXML::updateDataContainerWeights(
    DataContainer &DataCon, const std::vector<ComPWA::Event> &DataSample) {
  if (!DataCon.Weights && DataCon.WeightSum == 0.0) {
    LOG(INFO) << "Setting phase space sample weights...";
    std::vector<double> DataSetWeights;
    DataSetWeights.reserve(DataSample.size());
    double WeightSum(0.0);
    bool UniformWeights(true);

    for (auto const &x : DataSample) {
      DataSetWeights.push_back(x.Weight);
      WeightSum += x.Weight;
      if (x.Weight != 1.0)
        UniformWeights = false;
    }
    if (!UniformWeights) {
      DataCon.Weights = FunctionTree::MDouble("Weight", DataSetWeights);
    }
    DataCon.WeightSum = WeightSum;
  }
}

void IntensityBuilderXML::addFunctionTreeComponent(
    std::string Name, std::string Type,
    std::shared_ptr<ComPWA::FunctionTree::FunctionTree> FT) {
  if (nullptr != FT) {
    auto InsertResult = UniqueComponentFTMapping.insert(
        std::make_pair(Name, std::make_pair(Type, FT)));
    if (!InsertResult.second) {
      LOG(ERROR) << "IntensityBuilderXML::addFunctionTreeComponent(): "
                    "FunctionTree with name "
                 << Name << " already exists!";
    }
  }
}

FourMomentum createFourMomentum(const boost::property_tree::ptree &pt) {
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

ParticleStateTransitionKinematicsInfo
createKinematicsInfo(const ComPWA::ParticleList &PartList,
                     const boost::property_tree::ptree &pt) {
  auto initialS = pt.get_child("InitialState");
  auto InitialState = std::vector<int>(initialS.size());
  unsigned int counter(0);
  for (auto i : initialS) {
    std::string name = i.second.get<std::string>("<xmlattr>.Name");
    auto partP = ComPWA::findParticle(PartList, name);
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
    auto partP = ComPWA::findParticle(PartList, name);
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
        InitialState, FinalState, PartList, InitialStateP4,
        FinalStateEventPositionMapping);
  }

  return ParticleStateTransitionKinematicsInfo(
      InitialState, FinalState, PartList, FinalStateEventPositionMapping);
}

HelicityKinematics
createHelicityKinematics(const ComPWA::ParticleList &PartList,
                         const boost::property_tree::ptree &pt) {
  auto kininfo = createKinematicsInfo(PartList, pt);

  auto phspVal = pt.get_optional<double>("PhspVolume");
  if (phspVal) {
    return HelicityKinematics(kininfo, phspVal.get());
  } else {
    return HelicityKinematics(kininfo);
  }
}

} // namespace Physics
} // namespace ComPWA
