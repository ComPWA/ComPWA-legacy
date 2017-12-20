// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <sstream>

#include "Physics/DecayDynamics/RelativisticBreitWigner.hpp"
#include "Physics/DecayDynamics/AmpFlatteRes.hpp"
#include "Physics/DecayDynamics/NonResonant.hpp"

#include "Physics/HelicityFormalism/HelicityDecay.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {
HelicityDecay::HelicityDecay(std::shared_ptr<PartList> partL,
                             std::shared_ptr<Kinematics> kin,
                             const boost::property_tree::ptree &pt)
    : SubSys(pt) {
  load(partL, kin, pt);
}

void HelicityDecay::load(std::shared_ptr<PartList> partL,
                         std::shared_ptr<Kinematics> kin,
                         const boost::property_tree::ptree &pt) {

  LOG(trace) << "HelicityDecay::load() |";
  SubSys = SubSystem(pt);
  DataPosition =
      3 * std::dynamic_pointer_cast<HelicityKinematics>(kin)->dataID(SubSys);
  setPhspVolume(kin->phspVolume());

 Name = pt.get<std::string>("<xmlattr>.Name", "empty");
  Magnitude =
      std::make_shared<ComPWA::FitParameter>("Magnitude_" +Name, 1.0);
  Phase = std::make_shared<ComPWA::FitParameter>("Phase_" +Name, 0.0);
  std::shared_ptr<FitParameter> mag, phase;
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter") {
      if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude") {
        Magnitude = std::make_shared<FitParameter>(v.second);
      }
      if (v.second.get<std::string>("<xmlattr>.Type") == "Phase") {
        Phase = std::make_shared<FitParameter>(v.second);
      }
    } else {
      // ignored further settings. Should we throw an error?
    }
  }

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  auto partItr = partL->find(name);
  if (partItr == partL->end())
    throw std::runtime_error("HelicityDecay::load | Particle " + name +
                             " not found in list!");
  ComPWA::Spin J = partItr->second.GetSpinQuantumNumber("Spin");
  ComPWA::Spin mu(pt.get<double>("DecayParticle.<xmlattr>.Helicity"));

  // Read name and helicities from decay products
  auto decayProducts = pt.get_child("DecayProducts");
  if (decayProducts.size() != 2)
    throw boost::property_tree::ptree_error(
        "HelicityDecay::load() | Expect exactly two decay products (" +
        std::to_string(decayProducts.size()) + " given)!");

  auto p = decayProducts.begin();
  DecayProducts.first = p->second.get<std::string>("<xmlattr>.Name");
  DecayHelicities.first =
      ComPWA::Spin(p->second.get<int>("<xmlattr>.Helicity"));
  ++p;
  DecayProducts.second = p->second.get<std::string>("<xmlattr>.Name");
  DecayHelicities.second =
      ComPWA::Spin(p->second.get<double>("<xmlattr>.Helicity"));

  // Two-body decay
  if (SubSys.GetFinalStates().size() == 2) {
    // Create WignerD object
    AngularDist = std::make_shared<HelicityFormalism::AmpWignerD>(
        J, mu, DecayHelicities.first - DecayHelicities.second);

    auto partProp = partL->find(name)->second;
    std::string decayType = partProp.GetDecayType();

    if (decayType == "stable") {
      throw std::runtime_error("HelicityDecay::Factory() | Stable particle is "
                               "given as mother particle of a decay. Makes no "
                               "sense!");
    } else if (decayType == "relativisticBreitWigner") {
      DynamicFcn = std::make_shared<DecayDynamics::RelativisticBreitWigner>(
          name, DecayProducts, partL);
    } else if (decayType == "flatte") {
      DynamicFcn = std::make_shared<DecayDynamics::AmpFlatteRes>(
          name, DecayProducts, partL);
    } else {
      throw std::runtime_error(
          "HelicityDecay::Factory() | Unknown decay type " + decayType + "!");
    }

    // make sure dynamical function is created and set first
  } else { // Multi-body decay
    DynamicFcn = std::make_shared<DecayDynamics::NonResonant>(name);
    // We assume the we have a multi-body decay and assume that the decay
    // proceeds via constant (non-resonant) dynamics
    AngularDist = std::make_shared<AmpWignerD>(ComPWA::Spin(0), ComPWA::Spin(0),
                                               ComPWA::Spin(0));
  }
}

boost::property_tree::ptree HelicityDecay::save() const {
  auto pt = SubSys.save();
  pt.put<std::string>("<xmlattr>.Name", name());

  boost::property_tree::ptree tmp = Magnitude->save();
  tmp.put("<xmlattr>.Type", "Magnitude");
  pt.add_child("Parameter", tmp);

  tmp = Phase->save();
  tmp.put("<xmlattr>.Type", "Phase");
  pt.add_child("Parameter", tmp);

  pt.put("DecayParticle.<xmlattr>.Name", DynamicFcn->name());
  pt.put("DecayParticle.<xmlattr>.Helicity", AngularDist->mu());

  // TODO: put helicities of daughter particles
  return pt;
}

bool HelicityDecay::isModified() const {
  if (PartialAmplitude::isModified())
    return true;
  if (DynamicFcn->isModified()) {
    const_cast<double &>(CurrentIntegral) = integral();
    DynamicFcn->setModified(false);
    return true;
  }
  return false;
}

double HelicityDecay::normalization() const {
  if (DynamicFcn->isModified() || !CurrentIntegral)
    const_cast<double &>(CurrentIntegral) = integral();
  DynamicFcn->setModified(false);
  assert(CurrentIntegral != 0.0);
  return 1 / std::sqrt(CurrentIntegral);
}

std::shared_ptr<FunctionTree>
HelicityDecay::tree(std::shared_ptr<Kinematics> kin,
                    const ParameterList &sample, const ParameterList &toySample,
                    std::string suffix) {

  size_t n = sample.mDoubleValue(0)->values().size();
  size_t phspSize = toySample.mDoubleValue(0)->values().size();

  std::string nodeName = "PartialAmplitude(" + name() + ")" + suffix;

  auto tr = std::make_shared<FunctionTree>(
      nodeName, MComplex("", n), std::make_shared<MultAll>(ParType::MCOMPLEX));
  tr->createNode("Strength", std::make_shared<Value<std::complex<double>>>(),
                 std::make_shared<Complexify>(ParType::COMPLEX), nodeName);
  tr->createLeaf("Magnitude", Magnitude, "Strength");
  tr->createLeaf("Phase", Phase, "Strength");
  tr->createLeaf("PreFactor", PreFactor, nodeName);
  tr->insertTree(AngularDist->tree(sample, DataPosition + 1, DataPosition + 2),
                 nodeName);
  tr->insertTree(DynamicFcn->tree(sample, DataPosition), nodeName);

  tr->createNode("Normalization", std::make_shared<Value<double>>(),
                 std::make_shared<Inverse>(ParType::DOUBLE),
                 nodeName); // 1/normLH
  tr->createNode("SqrtIntegral", std::make_shared<Value<double>>(),
                 std::make_shared<SquareRoot>(ParType::DOUBLE),
                 "Normalization");
  tr->createNode("Integral", std::make_shared<Value<double>>(),
                 std::make_shared<MultAll>(ParType::DOUBLE), "SqrtIntegral");
  tr->createLeaf("PhspVolume", PhspVolume, "Integral");
  tr->createLeaf("InverseSampleSize", 1 / ((double)phspSize), "Integral");
  tr->createNode("Sum", std::make_shared<Value<double>>(),
                 std::make_shared<AddAll>(ParType::DOUBLE), "Integral");
  tr->createNode("Intensity", MDouble("", phspSize),
                 std::make_shared<AbsSquare>(ParType::MDOUBLE),
                 "Sum"); //|T_{ev}|^2
  tr->createNode("mult", MComplex("", phspSize),
                 std::make_shared<MultAll>(ParType::MCOMPLEX), "Intensity");
  tr->insertTree(
      AngularDist->tree(toySample, DataPosition + 1, DataPosition + 2, "_norm"),
      "mult" + suffix);
  tr->insertTree(DynamicFcn->tree(toySample, DataPosition, "_norm"),
                 "mult" + suffix);

  tr->parameter();
  return tr;
}

void HelicityDecay::parameters(ParameterList &list) {
  PartialAmplitude::parameters(list);
  DynamicFcn->parameters(list);
}

void HelicityDecay::updateParameters(const ParameterList &list) {
  PartialAmplitude::updateParameters(list);
  DynamicFcn->updateParameters(list);

  return;
}

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA
