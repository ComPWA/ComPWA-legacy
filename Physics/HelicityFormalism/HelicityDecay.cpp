// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <sstream>

#include "Physics/DecayDynamics/AmpFlatteRes.hpp"
#include "Physics/DecayDynamics/NonResonant.hpp"
#include "Physics/DecayDynamics/RelativisticBreitWigner.hpp"
#include "Physics/DecayDynamics/Voigtian.hpp"
#include "Physics/qft++/WignerD.h"

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

  LOG(TRACE) << "HelicityDecay::load() |";
  SubSys = SubSystem(pt);
  DataPosition =
      3 *
      std::dynamic_pointer_cast<HelicityKinematics>(kin)->addSubSystem(SubSys);
  setPhspVolume(kin->phspVolume());

  Name = pt.get<std::string>("<xmlattr>.Name", "empty");
  Magnitude = std::make_shared<ComPWA::FitParameter>("Magnitude_" + Name, 1.0);
  Phase = std::make_shared<ComPWA::FitParameter>("Phase_" + Name, 0.0);
  std::shared_ptr<FitParameter> mag, phase;
  PreFactor = std::complex<double>(1, 0);
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter") {
      if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude") {
        Magnitude = std::make_shared<FitParameter>(v.second);
      }
      if (v.second.get<std::string>("<xmlattr>.Type") == "Phase") {
        Phase = std::make_shared<FitParameter>(v.second);
      }
    } else if (v.first == "PreFactor") {
      double r = v.second.get<double>("<xmlattr>.Magnitude");
      double p = v.second.get<double>("<xmlattr>.Phase");
      PreFactor = std::polar(r, p);
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
  // if the node OrbitalAngularMomentum does not exist, set it to spin J as
  // default value
  ComPWA::Spin orbitL(pt.get<double>(
      "DecayParticle.<xmlattr>.OrbitalAngularMomentum", (double)J));
  
  const auto &canoSum = pt.get_child_optional("CanonicalSum");
  if (canoSum) {
    const auto &sumTree = canoSum.get();
    orbitL = sumTree.get<double>("<xmlattr>.L");
    double coef = sqrt((2 * (double) orbitL + 1) / (2 * (double) J + 1)) ;
    for (const auto &cg : sumTree.get_child("")) {
      if (cg.first != "ClebschGorden") continue;
      double j1 = cg.second.get<double>("<xmlattr>.j1");
      double m1 = cg.second.get<double>("<xmlattr>.m1");
      double j2 = cg.second.get<double>("<xmlattr>.j2");
      double m2 = cg.second.get<double>("<xmlattr>.m2");
      double J = cg.second.get<double>("<xmlattr>.J");
      double M = cg.second.get<double>("<xmlattr>.M");
      coef *= ComPWA::Physics::QFT::Clebsch(j1, m1, j2, m2, J, M);
    }
    PreFactor *= coef;
  }

  // Read name and helicities from decay products
  auto decayProducts = pt.get_child("DecayProducts");
  if (decayProducts.size() != 2)
    throw boost::property_tree::ptree_error(
        "HelicityDecay::load() | Expect exactly two decay products (" +
        std::to_string(decayProducts.size()) + " given)!");

  auto p = decayProducts.begin();
  DecayProducts.first = p->second.get<std::string>("<xmlattr>.Name");
  DecayHelicities.first =
      ComPWA::Spin(p->second.get<double>("<xmlattr>.Helicity"));
  ++p;
  DecayProducts.second = p->second.get<std::string>("<xmlattr>.Name");
  DecayHelicities.second =
      ComPWA::Spin(p->second.get<double>("<xmlattr>.Helicity"));

  // Two-body decay
  if (SubSys.getFinalStates().size() == 2) {
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
      DynamicFcn->SetOrbitalAngularMomentum(orbitL);
    } else if (decayType == "flatte") {
      DynamicFcn = std::make_shared<DecayDynamics::AmpFlatteRes>(
          name, DecayProducts, partL);
      DynamicFcn->SetOrbitalAngularMomentum(orbitL);
    } else if (decayType == "voigt") {
      DynamicFcn =
          std::make_shared<DecayDynamics::Voigtian>(name, DecayProducts, partL);
    } else if (decayType == "virtual" || decayType == "nonResonant") {
      DynamicFcn = std::make_shared<DecayDynamics::NonResonant>(name);
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
  pt.put("DecayParticle.<xmlatrr>.OrbitalAngularMomentum",
         DynamicFcn->GetOrbitalAngularMomentum());

  // TODO: put helicities of daughter particles
  return pt;
}

bool HelicityDecay::isModified() const {
  if (DynamicFcn->isModified() || magnitude() != CurrentMagnitude ||
      phase() != CurrentPhase) {
    return true;
  }
  return false;
}

void HelicityDecay::setModified(bool b) {
  DynamicFcn->setModified(b);
  if (b) {
    const_cast<double &>(CurrentIntegral) =
        std::numeric_limits<double>::quiet_NaN();
    const_cast<double &>(CurrentMagnitude) =
        std::numeric_limits<double>::quiet_NaN();
    const_cast<double &>(CurrentPhase) =
        std::numeric_limits<double>::quiet_NaN();
  } else {
    const_cast<double &>(CurrentIntegral) = integral();
    const_cast<double &>(CurrentMagnitude) = magnitude();
    const_cast<double &>(CurrentPhase) = phase();
  }
}

double HelicityDecay::normalization() const {
  if (isModified()) {
    // We have to do an ugly const cast here. Otherwise we would have to remove
    // all constness from all memeber functions...The basic design problem here
    // is that member variables are (smart) pointer which can be changed from
    // outside.
    const_cast<HelicityDecay *>(this)->setModified(false);
  }

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
