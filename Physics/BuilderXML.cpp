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
#include "Physics/Dynamics/RelativisticBreitWigner.hpp"
#include "Physics/Dynamics/Voigtian.hpp"
#include "Physics/HelicityFormalism/WignerD.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

namespace ComPWA {
namespace Physics {

using ComPWA::FunctionTree::FitParameter;
using ComPWA::FunctionTree::FunctionTreeIntensity;

IntensityBuilderXML::IntensityBuilderXML(
    ParticleList ParticleList, Kinematics &Kinematics,
    const boost::property_tree::ptree &ModelTree,
    const EventCollection &TruePhspSample,
    const EventCollection &RecoPhspSample)
    : PartList(ParticleList), Kinematic(Kinematics), ModelTree(ModelTree),
      TruePhspSample(TruePhspSample),
      RecoPhspSample(
          [&TruePhspSample, &RecoPhspSample]() -> const EventCollection & {
            if (TruePhspSample.Events.size() > 0 &&
                RecoPhspSample.Events.size() == 0)
              return TruePhspSample;
            else {
              return RecoPhspSample;
            }
          }()) {}

ComPWA::FunctionTree::FunctionTreeIntensity
IntensityBuilderXML::createIntensity() {
  LOG(TRACE) << "loading intensity...";
  // BlankState
  Parameters = FunctionTree::ParameterList();
  PhspData = DataContainer();
  PhspRecoData = DataContainer();

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> FT =
      createIntensityFT(ModelTree, Data.Data);

  updateDataContainerContent();

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
    std::map<
        std::string,
        std::pair<std::string, std::shared_ptr<ComPWA::FunctionTree::TreeNode>>>
        NewUniqueComponentFTMapping;
    std::string Type("");
    std::vector<std::shared_ptr<ComPWA::FunctionTree::TreeNode>> FTList;
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
              createCoherentIntensityFT(FTList), Parameters, Data.Data)));
      LOG(INFO) << "as a CoherentIntensity";
    } else {
      IntensityComponents.push_back(std::make_pair(
          ComponentName,
          std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensity>(
              createIncoherentIntensityFT(FTList), Parameters, Data.Data)));
      LOG(INFO) << "as a IncoherentIntensity";
    }
  }
  return IntensityComponents;
}

std::map<std::string, std::string>
IntensityBuilderXML::getAllComponentNames() const {
  std::map<std::string, std::string> Names;
  for (auto x : UniqueComponentFTMapping) {
    Names[x.first] = x.second.first;
  }
  return Names;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createIntensityFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample) {
  LOG(TRACE) << "loading intensity...";

  std::string IntensityClass(pt.get<std::string>("<xmlattr>.Class"));

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> FT(nullptr);

  if (IntensityClass == "IncoherentIntensity") {
    FT = createIncoherentIntensityFT(pt, DataSample);
  } else if (IntensityClass == "CoherentIntensity") {
    FT = createCoherentIntensityFT(pt, DataSample);
  } else if (IntensityClass == "StrengthIntensity") {
    FT = createStrengthIntensityFT(pt, DataSample);
  } else if (IntensityClass == "NormalizedIntensity") {
    FT = createNormalizedIntensityFT(pt, DataSample);
  } else {
    throw BadConfig("IntensityBuilderXML::createIntensityFT() | Found "
                    "unknown intensity " +
                    IntensityClass);
  }

  auto Component = pt.get_optional<std::string>("<xmlattr>.Component");
  if (Component.is_initialized()) {
    addFunctionTreeComponent(Component.get(), "Intensity", FT);
  }

  return FT;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createIncoherentIntensityFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample) {
  LOG(TRACE) << "constructing IncoherentIntensity ...";

  std::vector<std::shared_ptr<ComPWA::FunctionTree::TreeNode>> intens;
  for (const auto &x : pt) {
    if (x.first == "Intensity") {
      intens.push_back(createIntensityFT(x.second, DataSample));
    }
  }

  return createIncoherentIntensityFT(intens);
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createIncoherentIntensityFT(
    std::vector<std::shared_ptr<ComPWA::FunctionTree::TreeNode>> Intensities) {
  LOG(TRACE) << "constructing IncoherentIntensity ...";

  using namespace ComPWA::FunctionTree;

  auto tr =
      std::make_shared<TreeNode>(std::make_shared<AddAll>(ParType::MDOUBLE));

  for (auto x : Intensities) {
    tr->addNode(x);
  }

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createCoherentIntensityFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample) {
  LOG(TRACE) << "constructing CoherentIntensity ...";

  std::vector<std::shared_ptr<ComPWA::FunctionTree::TreeNode>> amps;
  for (const auto &x : pt) {
    if (x.first == "Amplitude") {
      amps.push_back(createAmplitudeFT(x.second, DataSample));
    }
  }

  return createCoherentIntensityFT(amps);
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createCoherentIntensityFT(
    std::vector<std::shared_ptr<ComPWA::FunctionTree::TreeNode>> Amplitudes) {
  LOG(TRACE) << "constructing CoherentIntensity ...";

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<ComPWA::FunctionTree::TreeNode>(
      std::make_shared<AbsSquare>(ParType::MDOUBLE));

  auto SumOfAmps =
      std::make_shared<TreeNode>(std::make_shared<AddAll>(ParType::MCOMPLEX));
  tr->addNode(SumOfAmps);
  for (auto x : Amplitudes) {
    SumOfAmps->addNode(x);
  }

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createStrengthIntensityFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample) {
  LOG(TRACE) << "creating StrengthIntensity ...";

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

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<ComPWA::FunctionTree::TreeNode>(
      MDouble("", 0), std::make_shared<MultAll>(ParType::MDOUBLE));

  tr->addNode(FunctionTree::createLeaf(Strength));
  tr->addNode(createIntensityFT(UndecoratedIntensityPT, DataSample));

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createNormalizedIntensityFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample) {
  LOG(TRACE) << "creating NormalizedIntensity ...";

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
        LOG(TRACE)
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

  return normalizeIntensityFT(UnnormalizedPT, DataSample, IntegratorClassName);
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::normalizeIntensityFT(
    const boost::property_tree::ptree &UnnormalizedPT,
    const ComPWA::FunctionTree::ParameterList &DataSample,
    std::string IntegratorClassName) {
  LOG(TRACE) << "creating Normalized FunctionTree ...";

  if (RecoPhspSample.Events.size() == 0)
    LOG(FATAL) << "IntensityBuilderXML::normalizeIntensityFT(): "
                  "reco phsp sample is not set!";
  updateDataContainerWeights(PhspRecoData, RecoPhspSample);

  using namespace ComPWA::FunctionTree;

  auto NormalizedFT = std::make_shared<ComPWA::FunctionTree::TreeNode>(
      std::make_shared<MultAll>(ParType::MDOUBLE));

  auto FTData = createIntensityFT(UnnormalizedPT, DataSample);
  NormalizedFT->addNode(FTData);

  bool PreviousSetting = ComponentRegisteringEnabled;
  ComponentRegisteringEnabled = false;
  auto FTPhspData = createIntensityFT(UnnormalizedPT, PhspRecoData.Data);
  ComponentRegisteringEnabled = PreviousSetting;
  auto normtree =
      createIntegrationStrategyFT(FTPhspData, PhspRecoData.Weights,
                                  PhspRecoData.WeightSum, IntegratorClassName);

  NormalizedFT->addNode(normtree);

  return NormalizedFT;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createIntegrationStrategyFT(
    std::shared_ptr<ComPWA::FunctionTree::TreeNode> UnnormalizedIntensity,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        PhspWeights,
    double PhspWeightSum, std::string IntegratorClassName) {
  LOG(TRACE) << "creating IntegrationStrategy ...";

  using namespace ComPWA::FunctionTree;
  std::shared_ptr<TreeNode> tr;

  if (IntegratorClassName == "MCIntegrationStrategy") {
    // update the PhspData container
    updateDataContainerState();

    tr = std::make_shared<TreeNode>(
        std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)));
    auto Integral = std::make_shared<TreeNode>(
        std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)));
    tr->addNode(Integral);
    Integral->addNodes(
        {createLeaf(Kinematic.phspVolume(), "PhspVolume"),
         createLeaf(1.0 / PhspWeightSum, "Inverse_PhspWeightSum")});
    auto Sum = std::make_shared<TreeNode>(
        std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)));
    Integral->addNode(Sum);
    auto WeightedIntensities = std::make_shared<TreeNode>(
        std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
    Sum->addNode(WeightedIntensities);

    if (PhspWeights)
      WeightedIntensities->addNode(createLeaf(PhspWeights));
    WeightedIntensities->addNode(UnnormalizedIntensity);
  } else {
    LOG(WARNING) << "IntensityBuilderXML::createIntegrationStrategyFT(): "
                    "IntegrationStrategy type "
                 << IntegratorClassName << " unknown!";
  }

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createAmplitudeFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample) {
  auto ampclass = pt.get<std::string>("<xmlattr>.Class");

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> FT(nullptr);

  if (ampclass == "HelicityDecay") {
    FT = createHelicityDecayFT(pt, DataSample);
  } else if (ampclass == "CoefficientAmplitude") {
    FT = createCoefficientAmplitudeFT(pt, DataSample);
  } else if (ampclass == "SequentialAmplitude") {
    FT = createSequentialAmplitudeFT(pt, DataSample);
  } else if (ampclass == "NormalizedAmplitude") {
    FT = createNormalizedAmplitudeFT(pt, DataSample);
  } else {
    throw BadConfig(
        "IntensityBuilderXML::createAmplitude(): Unknown amplitude " +
        ampclass);
  }

  auto Component = pt.get_optional<std::string>("<xmlattr>.Component");
  if (Component.is_initialized()) {
    addFunctionTreeComponent(Component.get(), "Amplitude", FT);
  }

  return FT;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createNormalizedAmplitudeFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample) {
  LOG(TRACE) << "creating NormalizedAmplitude ...";

  if (TruePhspSample.Events.size() == 0)
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
        LOG(TRACE)
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

  auto FTData = createAmplitudeFT(UnnormalizedPT, DataSample);

  using namespace ComPWA::FunctionTree;

  auto NormalizedFT =
      std::make_shared<TreeNode>(std::make_shared<MultAll>(ParType::MCOMPLEX));

  NormalizedFT->addNode(FTData);

  bool PreviousSetting = ComponentRegisteringEnabled;
  ComponentRegisteringEnabled = false;
  auto FTPhspData = createAmplitudeFT(UnnormalizedPT, PhspData.Data);
  ComponentRegisteringEnabled = PreviousSetting;
  // this phspdata function tree has to be made into a double valued function
  auto FTPhspDataAbsSquared =
      std::make_shared<TreeNode>(std::make_shared<AbsSquare>(ParType::MDOUBLE));
  FTPhspDataAbsSquared->addNode(FTPhspData);

  auto normtreesquared =
      createIntegrationStrategyFT(FTPhspDataAbsSquared, PhspData.Weights,
                                  PhspData.WeightSum, IntegratorClassName);

  auto normtree =
      std::make_shared<TreeNode>(std::make_shared<SquareRoot>(ParType::DOUBLE));

  normtree->addNode(normtreesquared);

  NormalizedFT->addNode(normtree);

  return NormalizedFT;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createCoefficientAmplitudeFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample) {
  LOG(TRACE) << "constructing CoefficientAmplitude ...";

  std::shared_ptr<FitParameter> Magnitude(nullptr);
  std::shared_ptr<FitParameter> Phase(nullptr);
  auto PreFactor = std::complex<double>(1, 0);
  bool IsPreFactorSet(false);
  boost::property_tree::ptree Amplitude;
  for (const auto &v : pt) {
    if (v.first == "Parameter") {
      if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude")
        Magnitude = std::make_shared<FitParameter>(v.second);
      if (v.second.get<std::string>("<xmlattr>.Type") == "Phase")
        Phase = std::make_shared<FitParameter>(v.second);
    } else if (v.first == "PreFactor") {
      IsPreFactorSet = true;
      boost::optional<double> optr =
          v.second.get_optional<double>("<xmlattr>.Magnitude");
      if (optr.is_initialized()) {
        double r(optr.value());
        if (r < 0.0)
          throw BadConfig(
              "IntensityBuilderXML::createCoefficientAmplitudeFT() | "
              "PreFactor Magnitude below zero!");
        double Phase(0.0);
        boost::optional<double> OptPhase =
            v.second.get_optional<double>("<xmlattr>.Phase");
        if (OptPhase.is_initialized())
          Phase = OptPhase.value();
        PreFactor = std::polar(r, Phase);
      } else {
        double Real = v.second.get<double>("<xmlattr>.Real");
        double Imaginary(0.0);
        boost::optional<double> OptIm =
            v.second.get_optional<double>("<xmlattr>.Imaginary");
        if (OptIm.is_initialized())
          Imaginary = OptIm.value();
        PreFactor = std::complex<double>(Real, Imaginary);
      }
    } else if (v.first == "Amplitude") {
      Amplitude = v.second;
    } else if (v.first != "<xmlattr>") {
      throw BadConfig("IntensityBuilderXML::createCoefficientAmplitudeFT() | "
                      "Unknown tag " +
                      v.first + "!");
    }
  }

  auto amp_ft = createAmplitudeFT(Amplitude, DataSample);

  auto tr = std::make_shared<ComPWA::FunctionTree::TreeNode>(
      ComPWA::FunctionTree::MComplex("", 0),
      std::make_shared<ComPWA::FunctionTree::MultAll>(
          ComPWA::FunctionTree::ParType::MCOMPLEX));

  bool AddStrength(true);
  if (IsPreFactorSet) {
    tr->addNode(FunctionTree::createLeaf(PreFactor, "PreFactor"));
    if (!Magnitude && !Phase) {
      AddStrength = false;
    } else if (!Magnitude && Phase) {
      throw BadParameter("IntensityBuilderXML::createCoefficientAmplitude() | "
                         "No magnitude parameter found while phase is set.");
    } else if (!Phase && Magnitude)
      throw BadParameter("IntensityBuilderXML::createCoefficientAmplitude() | "
                         "No phase parameter found, while magnitude is set.");
  } else {
    if (!Magnitude)
      throw BadParameter("IntensityBuilderXML::createCoefficientAmplitude() | "
                         "No magnitude parameter found and no prefactor set.");
    if (!Phase)
      throw BadParameter("IntensityBuilderXML::createCoefficientAmplitude() | "
                         "No phase parameter found, while magnitude is set.");
  }

  if (AddStrength) {
    Magnitude = Parameters.addUniqueParameter(Magnitude);
    Phase = Parameters.addUniqueParameter(Phase);

    auto Strength = std::make_shared<ComPWA::FunctionTree::TreeNode>(
        std::make_shared<ComPWA::FunctionTree::Value<std::complex<double>>>(),
        std::make_shared<ComPWA::FunctionTree::Complexify>(
            ComPWA::FunctionTree::ParType::COMPLEX));
    Strength->addNodes({createLeaf(Magnitude), createLeaf(Phase)});

    tr->addNode(Strength);
  }
  tr->addNode(amp_ft);

  return tr;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createSequentialAmplitudeFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample) {
  LOG(TRACE) << "constructing SequentialAmplitude ...";

  std::vector<std::shared_ptr<ComPWA::FunctionTree::TreeNode>> Amplitudes;
  for (const auto &v : pt) {
    if (v.first == "Amplitude") {
      std::shared_ptr<ComPWA::FunctionTree::TreeNode> AmpTree =
          createAmplitudeFT(v.second, DataSample);

      Amplitudes.push_back(AmpTree);
    } else if (v.first != "<xmlattr>") {
      throw BadConfig("SequentialAmplitude::createSequentialAmplitude() | "
                      "Unknown tag " +
                      v.first + "!");
    }
  }

  using namespace ComPWA::FunctionTree;
  auto tr =
      std::make_shared<TreeNode>(std::make_shared<MultAll>(ParType::MCOMPLEX));

  for (auto x : Amplitudes) {
    tr->addNode(x);
  }

  return tr;
}

TwoBodyDecayInfo extractDecayInfo(const boost::property_tree::ptree &pt) {
  // Read name and helicities from decay products
  auto decayProducts = pt.get_child("DecayProducts");
  if (decayProducts.size() != 2)
    throw boost::property_tree::ptree_error(
        "IntensityBuilderXML::createHelicityDecayFT(): Expect exactly two "
        "decay products (" +
        std::to_string(decayProducts.size()) + " given)!");

  auto recoil =
      pt.get_optional<std::string>("RecoilSystem.<xmlattr>.RecoilFinalState");
  IndexList RecoilState;
  if (recoil) {
    RecoilState = stringToVectInt(recoil.get());
  }
  auto parent_recoil = pt.get_optional<std::string>(
      "RecoilSystem.<xmlattr>.ParentRecoilFinalState");
  IndexList ParentRecoilState;
  if (parent_recoil) {
    ParentRecoilState = stringToVectInt(parent_recoil.get());
  }
  std::vector<IndexList> DecayProductStates;
  std::vector<std::string> Names;
  std::vector<double> Helicities;
  for (auto p : decayProducts) {
    DecayProductStates.push_back(
        stringToVectInt(p.second.get<std::string>("<xmlattr>.FinalState")));
    Names.push_back(p.second.get<std::string>("<xmlattr>.Name"));
    Helicities.push_back(p.second.get<double>("<xmlattr>.Helicity"));
  }
  return TwoBodyDecayInfo{
      SubSystem(DecayProductStates, RecoilState, ParentRecoilState),
      std::make_pair(Names.at(0), Names.at(1)),
      std::make_pair(Helicities.at(0), Helicities.at(1))};
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
IntensityBuilderXML::createHelicityDecayFT(
    const boost::property_tree::ptree &pt,
    const ComPWA::FunctionTree::ParameterList &DataSample) {
  LOG(TRACE) << "IntensityBuilderXML::createHelicityDecayFT(): ";
  auto &kin = dynamic_cast<HelicityFormalism::HelicityKinematics &>(Kinematic);
  auto DecayProductInfo = extractDecayInfo(pt);
  auto KinVarNames = kin.registerSubSystem(DecayProductInfo.SubSys);
  auto const &FinalStates = DecayProductInfo.SubSys.getFinalStates();
  std::pair<std::string, std::string> DecayProductInvMassNames;
  if (FinalStates.at(0).size() > 1) {
    DecayProductInvMassNames.first =
        kin.registerInvariantMassSquared(FinalStates.at(0));
  }
  if (FinalStates.at(1).size() > 1) {
    DecayProductInvMassNames.second =
        kin.registerInvariantMassSquared(FinalStates.at(1));
  }

  updateDataContainerState();

  auto InvMassSq =
      FunctionTree::findMDoubleValue(std::get<0>(KinVarNames), DataSample);
  auto Theta =
      FunctionTree::findMDoubleValue(std::get<1>(KinVarNames), DataSample);
  auto Phi =
      FunctionTree::findMDoubleValue(std::get<2>(KinVarNames), DataSample);

  auto parMass1 = std::make_shared<FitParameter>(
      ComPWA::findParticle(PartList, DecayProductInfo.Names.first).getMass());
  auto parMass2 = std::make_shared<FitParameter>(
      ComPWA::findParticle(PartList, DecayProductInfo.Names.second).getMass());
  parMass1 = Parameters.addUniqueParameter(parMass1);
  parMass2 = Parameters.addUniqueParameter(parMass2);

  auto InvMassDaughter1 = FunctionTree::MDouble(
      parMass1->name(), std::vector<double>{std::pow(parMass1->value(), 2)});
  if (DecayProductInvMassNames.first != "") {
    InvMassDaughter1 = FunctionTree::findMDoubleValue(
        DecayProductInvMassNames.first, DataSample);
  }
  auto InvMassDaughter2 = FunctionTree::MDouble(
      parMass2->name(), std::vector<double>{std::pow(parMass2->value(), 2)});
  if (DecayProductInvMassNames.second != "") {
    InvMassDaughter2 = FunctionTree::findMDoubleValue(
        DecayProductInvMassNames.second, DataSample);
  }

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  auto partProp = ComPWA::findParticle(PartList, name);

  double J = partProp.getQuantumNumber<double>("Spin");
  double mu(pt.get<double>("DecayParticle.<xmlattr>.Helicity"));
  // if the node OrbitalAngularMomentum does not exist, set it to spin J as
  // default value
  unsigned int L(J);

  auto Mass = std::make_shared<FunctionTree::FitParameter>(
      partProp.getMass().Name, partProp.getMass().Value,
      partProp.getMass().Error.first);
  Mass->fixParameter(partProp.getMass().IsFixed);
  Mass = Parameters.addUniqueParameter(Mass);

  double PreFactor(1.0);
  const auto &canoSum = pt.get_child_optional("CanonicalSum");
  if (canoSum) {
    const auto &sumTree = canoSum.get();
    L = sumTree.get<unsigned int>("<xmlattr>.L");
    double Coefficient = std::sqrt((2.0 * L + 1) / (2 * J + 1));
    for (const auto &cg : sumTree.get_child("")) {
      if (cg.first != "ClebschGordan")
        continue;
      double j1 = cg.second.get<double>("<xmlattr>.j1");
      double m1 = cg.second.get<double>("<xmlattr>.m1");
      double j2 = cg.second.get<double>("<xmlattr>.j2");
      double m2 = cg.second.get<double>("<xmlattr>.m2");
      double J = cg.second.get<double>("<xmlattr>.J");
      double M = cg.second.get<double>("<xmlattr>.M");
      Coefficient *= ComPWA::QFT::Clebsch(j1, m1, j2, m2, J, M);
    }
    PreFactor *= Coefficient;
  }

  auto decayInfo = partProp.getDecayInfo();
  std::string FFType("");
  std::shared_ptr<Dynamics::FormFactor> FormFactor =
      std::make_shared<Dynamics::NoFormFactor>();
  std::shared_ptr<FitParameter> parRadius;
  std::shared_ptr<FitParameter> Width;
  for (const auto &node : decayInfo.get_child("")) {
    if (node.first == "FormFactor") {
      FFType = node.second.get<std::string>("<xmlattr>.Type");
      if ("BlattWeisskopf" == FFType) {
        if (L > 4)
          throw std::runtime_error(
              "createHelicityDecayFT() | Blatt-Weisskopf form factors are "
              "implemented for L up to 4!");
        FormFactor = std::make_shared<Dynamics::BlattWeisskopfFormFactor>();
      } else if ("CrystalBarrel" == FFType) {
        if (L != 0)
          throw std::runtime_error(
              "createHelicityDecayFT() | Crystal Barrel form factors are "
              "implemented for L=0 only!");
        FormFactor = std::make_shared<Dynamics::CrystalBarrelFormFactor>();
      } else {
        // if user specifies unknown form factor
        throw std::runtime_error("createHelicityDecayFT() | Form factor type " +
                                 FFType + " not specified!");
      }
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

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> DynamicFunctionFT(nullptr);

  if (decayType == "stable") {
    throw std::runtime_error(
        "IntensityBuilderXML::createHelicityDecayFT(): Stable particle is "
        "given as mother particle of a decay. Makes no sense!");
  } else if (decayType.find("BreitWigner") != std::string::npos) {
    using namespace ComPWA::Physics::Dynamics::RelativisticBreitWigner;
    InputInfo RBW;
    RBW.Type = decayType;
    RBW.Mass = Mass;
    RBW.Width = Width;
    RBW.MesonRadius = parRadius;
    RBW.DaughterInvariantMasses =
        std::make_pair(InvMassDaughter1, InvMassDaughter2);
    RBW.FormFactorFunctor = FormFactor;
    RBW.L = L;
    DynamicFunctionFT = createFunctionTree(RBW, InvMassSq);
  } else if (decayType == "flatte") {
    ComPWA::Physics::Dynamics::Flatte::InputInfo FlatteInfo;
    FlatteInfo.Mass = Mass;
    FlatteInfo.MesonRadius = parRadius;
    FlatteInfo.DaughterInvariantMasses =
        std::make_pair(InvMassDaughter1, InvMassDaughter2);
    FlatteInfo.FormFactorFunctor = FormFactor;
    FlatteInfo.L = L;
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
    DynamicFunctionFT = createFunctionTree(FlatteInfo, InvMassSq);
  } else if (decayType == "voigt") {
    using namespace ComPWA::Physics::Dynamics::Voigtian;
    InputInfo VoigtInfo;
    VoigtInfo.Mass = Mass;
    VoigtInfo.Width = Width;
    VoigtInfo.MesonRadius = parRadius;
    VoigtInfo.DaughterInvariantMasses =
        std::make_pair(InvMassDaughter1, InvMassDaughter2);
    VoigtInfo.FormFactorFunctor = FormFactor;
    VoigtInfo.L = L;
    VoigtInfo.Sigma = decayInfo.get<double>("Resolution.<xmlattr>.Sigma");
    DynamicFunctionFT = createFunctionTree(VoigtInfo, InvMassSq);
  } else if (decayType == "virtual" || decayType == "nonResonant") {
    DynamicFunctionFT = FunctionTree::createLeaf(1, "Dynamics");
  } else {
    throw std::runtime_error("HelicityDecay::Factory() | Unknown decay type " +
                             decayType + "!");
  }

  auto AngularFunction =
      ComPWA::Physics::HelicityFormalism::WignerD::createFunctionTree(
          J, mu,
          DecayProductInfo.Helicities.first -
              DecayProductInfo.Helicities.second,
          Theta, Phi);

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<ComPWA::FunctionTree::TreeNode>(
      std::make_shared<MultAll>(ParType::MCOMPLEX));
  tr->addNodes(
      {createLeaf(PreFactor, "PreFactor"), AngularFunction, DynamicFunctionFT});

  // set production formfactor
  if (FFType != "" && L > 0) {
    if (parRadius == nullptr) {
      throw std::runtime_error("IntensityBuilderXML::createHelicityDecayFT(): "
                               "No MesonRadius given for amplitude " +
                               name +
                               "! It is needed to calculate the form factor!");
    }

    std::shared_ptr<ComPWA::FunctionTree::TreeNode> ProductionFormFactorFT =
        Dynamics::createFunctionTree(InvMassDaughter1, InvMassDaughter2,
                                     parRadius, L, FormFactor, InvMassSq);

    tr->addNode(ProductionFormFactorFT);
  }

  return tr;
}

void updateDataContainerState(ComPWA::FunctionTree::ParameterList &DataSample,
                              const Kinematics &Kinematics) {
  auto ExemplaryDataSet = Kinematics.convert({Kinematics.getFinalStatePIDs()});

  for (auto const &x : ExemplaryDataSet.Data) {
    std::vector<double> temp;
    // check if this key already exists than skip otherwise insert
    bool Exists(false);
    for (auto MDVal : DataSample.mDoubleValues()) {
      if (MDVal->name() == x.first) {
        Exists = true;
        break;
      }
    }
    if (!Exists) {
      DataSample.addValue(
          std::make_shared<ComPWA::FunctionTree::Value<std::vector<double>>>(
              x.first, temp));
    }
  }
}

void IntensityBuilderXML::updateDataContainerState() {
  LOG(TRACE) << "updating data containers...";
  Physics::updateDataContainerState(Data.Data, Kinematic);
  Physics::updateDataContainerState(PhspData.Data, Kinematic);
  Physics::updateDataContainerState(PhspRecoData.Data, Kinematic);
}

void updateDataContainerContent(ComPWA::FunctionTree::ParameterList &DataList,
                                const EventCollection &DataSample,
                                const Kinematics &Kinematics) {
  LOG(INFO) << "Updating data container content...";
  if (DataSample.Events.size() == 0)
    return;

  auto DataSet = Kinematics.convert(DataSample);

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
      DataList.mDoubleValue(i)->setValue(
          DataSet.Data[DataList.mDoubleValue(i)->name()]);
    }
  }
}

void IntensityBuilderXML::updateDataContainerContent() {
  Physics::updateDataContainerContent(PhspRecoData.Data, RecoPhspSample,
                                      Kinematic);
}

void IntensityBuilderXML::updateDataContainerWeights(
    DataContainer &DataCon, const EventCollection &DataSample) {
  if (!DataCon.Weights && DataCon.WeightSum == 0.0) {
    LOG(INFO) << "Setting phase space sample weights...";
    std::vector<double> DataSetWeights;
    DataSetWeights.reserve(DataSample.Events.size());
    double WeightSum(0.0);
    bool UniformWeights(true);

    for (auto const &Event : DataSample.Events) {
      DataSetWeights.push_back(Event.Weight);
      WeightSum += Event.Weight;
      if (Event.Weight != 1.0)
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
    std::shared_ptr<ComPWA::FunctionTree::TreeNode> FT) {
  if (nullptr != FT && ComponentRegisteringEnabled) {
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
  auto InitialState = std::vector<pid>(initialS.size());
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
  auto FinalState = std::vector<pid>(finalS.size());
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

HelicityKinematics createHelicityKinematics(const std::string XmlFile) {
  auto ParticleList = readParticles(XmlFile);
  return createHelicityKinematics(ParticleList, XmlFile);
}

HelicityKinematics
createHelicityKinematics(const ComPWA::ParticleList &ParticleList,
                         const std::string XmlFile) {
  boost::property_tree::ptree ptree;
  boost::property_tree::xml_parser::read_xml(XmlFile, ptree);
  auto it = ptree.find("HelicityKinematics");
  if (it != ptree.not_found()) {
    return ComPWA::Physics::createHelicityKinematics(ParticleList, it->second);
  } else {
    throw ComPWA::BadConfig("ComPWA::Physics::createHelicityKinematics(): "
                            "HelicityKinematics tag not found in xml file!");
  }
}

HelicityKinematics
createHelicityKinematics(const ComPWA::ParticleList &ParticleList,
                         const boost::property_tree::ptree &ptree) {
  auto kininfo = createKinematicsInfo(ParticleList, ptree);

  auto phspVal = ptree.get_optional<double>("PhspVolume");
  if (phspVal) {
    return HelicityKinematics(kininfo, phspVal.get());
  } else {
    return HelicityKinematics(kininfo);
  }
}

} // namespace Physics
} // namespace ComPWA
