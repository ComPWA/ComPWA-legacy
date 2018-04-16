// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <memory>
#include "Physics/HelicityFormalism/HelicityFromCanonicalSum.hpp"
#include <cmath>
#include "Physics/qft++/WignerD.h"

#include <string>

namespace ComPWA{
namespace Physics {
namespace HelicityFormalism {

HelicityFromCanonicalSum::HelicityFromCanonicalSum(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) {
  load(partL, kin, pt);
}

void HelicityFromCanonicalSum::load(std::shared_ptr<PartList> partL,
                                      std::shared_ptr<Kinematics> kin,
                                      const boost::property_tree::ptree &pt) {
  LOG(trace) << "HelicityFromCanonicalSum::Factory() | Construction....";
//  setName(pt.get<std::string>("<xmlattr>.Name", "empty"));

  boost::property_tree::ptree helDecTree(pt);
  helDecTree.erase("CanonicalSum");

  std::string helDecName = pt.get<std::string>("<xmlattr>.Name", "empty");
  std::string mothName = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  auto partItr = partL->find(mothName);
  if (partItr == partL->end())
    throw std::runtime_error("HelicityFromCanonicalSum::load | Particle " + mothName +
        " not found in list!");
  double spinJ = (double) partItr->second.GetSpinQuantumNumber("Spin");

  const auto &optSumTree = pt.get_child_optional("CanonicalSum");
  if (optSumTree) {
    auto const &sumTree = optSumTree.get();

    double orbitL = sumTree.get<double>("<xmlattr>.L");
    helDecTree.put("DecayParticle.<xmlattr>.OrbitalAngularMomentum",
        sumTree.get<std::string>("<xmlattr>.L"));
    helDecName += "_L";
    helDecName += sumTree.get<std::string>("<xmlattr>.L");
    helDecName += "_S";
    helDecName += sumTree.get<std::string>("<xmlattr>.S");
    helDecTree.put("<xmlattr>.Name", helDecName);

    std::shared_ptr<HelicityDecay> ampHelDec = std::make_shared<HelicityDecay>(partL, kin, helDecTree); 

    double normCoef = sqrt( (2 * orbitL + 1)/(2 * spinJ + 1) );
    double cgCoef = 1.0;
    for (const auto &daug : sumTree.get_child("")) {
      if (daug.first != "ClebschGorden") continue;
      double j1 = daug.second.get<double>("<xmlattr>.j1");
      double m1 = daug.second.get<double>("<xmlattr>.m1");
      double j2 = daug.second.get<double>("<xmlattr>.j2");
      double m2 = daug.second.get<double>("<xmlattr>.m2");
      double J = daug.second.get<double>("<xmlattr>.J");
      double M = daug.second.get<double>("<xmlattr>.M");
      cgCoef *= ComPWA::Physics::QFT::Clebsch(j1, m1, j2, m2, J, M);
    }
    PreFactor = std::complex<double>(normCoef * cgCoef, 0);
    ampHelDec->setPrefactor(PreFactor);

    addPartialAmplitude(ampHelDec);
  } else {
    LOG(trace) << "HelicityFromCanonicalSum::Factory() | <CanonicalSum /> in xml is not set!";
  }
}

boost::property_tree::ptree HelicityFromCanonicalSum::save() const {

  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", name());

  boost::property_tree::ptree tmp = Magnitude->save();
  tmp.put("<xmlattr>.Type", "Magnitude");
  pt.add_child("Parameter", tmp);

  tmp = Phase->save();
  tmp.put("<xmlattr>.Type", "Phase");
  pt.add_child("Parameter", tmp);

  auto pref = preFactor();
  if (pref != std::complex<double>(1, 0)) {
    boost::property_tree::ptree ppp;
    ppp.put("<xmlattr>.Magnitude", std::abs(pref));
    ppp.put("<xmlattr>.Phase", std::arg(pref));
    pt.add_child("PreFactor", ppp);
  }

  for (auto i : PartialAmplitudes) {
    pt.add_child("PartialAmplitude", i->save());
  }
  return pt;
}

std::shared_ptr<ComPWA::FunctionTree> HelicityFromCanonicalSum::tree(
    std::shared_ptr<Kinematics> kin, const ParameterList &sample,
    const ParameterList &toySample, std::string suffix) {

  size_t n = sample.mDoubleValue(0)->values().size();

  auto tr = std::make_shared<FunctionTree>(
      "Amplitude(" + name() + ")" + suffix, MComplex("", n),
      std::make_shared<MultAll>(ParType::MCOMPLEX));
  tr->createNode("Strength", std::make_shared<Value<std::complex<double>>>(),
                 std::make_shared<Complexify>(ParType::COMPLEX),
                 "Amplitude(" + name() + ")" + suffix);
  tr->createLeaf("Magnitude", Magnitude, "Strength");
  tr->createLeaf("Phase", Phase, "Strength");
  tr->createLeaf("PreFactor", preFactor(),
                 "Amplitude(" + name() + ")" + suffix);

  for (auto i : PartialAmplitudes) {
    std::shared_ptr<FunctionTree> resTree = i->tree(kin, sample, toySample, "");
    if (!resTree->sanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Amplitude tree didn't pass sanity check!");
    resTree->parameter();
    tr->insertTree(resTree, "Amplitude(" + name() + ")" + suffix);
  }

  return tr;
}

void HelicityFromCanonicalSum::parameters(ParameterList &list) {
  Amplitude::parameters(list);
  for (auto i : PartialAmplitudes) {
    i->parameters(list);
  }
}

void HelicityFromCanonicalSum::updateParameters(const ParameterList &list) {
  // Try to update magnitude
  std::shared_ptr<FitParameter> mag;
  try {
    mag = FindParameter(Magnitude->name(), list);
  } catch (std::exception &ex) {
  }
  if (mag)
    Magnitude->updateParameter(mag);
  std::shared_ptr<FitParameter> phase;

  // Try to update phase
  try {
    phase = FindParameter(Phase->name(), list);
  } catch (std::exception &ex) {
  }
  if (phase)
    Phase->updateParameter(phase);

  for (auto i : PartialAmplitudes)
    i->updateParameters(list);

  return;
}

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA
