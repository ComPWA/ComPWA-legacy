// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <memory>
#include "Physics/SequentialPartialAmplitude.hpp"
#include "Physics/HelicityFormalism/HelicityFromCanonicalSum.hpp"

using namespace ComPWA::Physics::HelicityFormalism;

SequentialPartialAmplitude::SequentialPartialAmplitude(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) {
  load(partL, kin, pt);
}

void SequentialPartialAmplitude::load(std::shared_ptr<PartList> partL,
                                      std::shared_ptr<Kinematics> kin,
                                      const boost::property_tree::ptree &pt) {
  LOG(trace) << "SequentialPartialAmplitude::Factory() | Construction....";
  setName(pt.get<std::string>("<xmlattr>.Name", "empty"));

  PreFactor = std::complex<double>(1, 0);
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter"){
      if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude")
        Magnitude = std::make_shared<FitParameter>(v.second);
      if (v.second.get<std::string>("<xmlattr>.Type") == "Phase")
        Phase = std::make_shared<FitParameter>(v.second);
    } else if (v.first == "PartialAmplitude" &&
               v.second.get<std::string>("<xmlattr>.Class") ==
                   "HelicityDecay") {
      addPartialAmplitude(
          std::make_shared<HelicityDecay>(partL, kin, v.second));
    } else if (v.first == "PartialAmplitude" &&
               v.second.get<std::string>("<xmlattr>.Class") ==
                   "HelicityFromCanonicalSum") {
      HelicityFromCanonicalSum canonicalSum(partL, kin, v.second);
      for (auto &i : canonicalSum.partialAmplitudes()) {
        addPartialAmplitude(i);
      }
    } else if (v.first == "PartialAmplitude" &&
               v.second.get<std::string>("<xmlattr>.Class") ==
                   "NonResonant") {
      addPartialAmplitude(
          std::make_shared<ComPWA::Physics::NonResonant>(partL, kin, v.second));
    } else if (v.first == "PreFactor") {
      double r = v.second.get<double>("<xmlattr>.Magnitude");
      double p = v.second.get<double>("<xmlattr>.Phase");
      PreFactor = std::polar(r, p);
    } else {
      // ignored further settings. Should we throw an error?
    }
  }

  if (!Magnitude)
    throw BadParameter("SequentialPartialAmplitude::Factory() | No magnitude "
                       "parameter found.");
  if (!Phase)
    throw BadParameter(
        "SequentialPartialAmplitude::Factory() | No phase parameter found.");
}

boost::property_tree::ptree SequentialPartialAmplitude::save() const {

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

bool SequentialPartialAmplitude::isModified() const {
  if (magnitude() != CurrentMagnitude || phase() != CurrentPhase)
    return true;

  for (auto i : PartialAmplitudes)
    if (i->isModified())
      return true;
  return false;
}

void SequentialPartialAmplitude::setModified(bool b) {

  if (b) {
    const_cast<double &>(CurrentMagnitude) =
        std::numeric_limits<double>::quiet_NaN();
    const_cast<double &>(CurrentPhase) =
        std::numeric_limits<double>::quiet_NaN();
  } else {
    const_cast<double &>(CurrentMagnitude) = magnitude();
    const_cast<double &>(CurrentPhase) = phase();
  }
  for (auto i : PartialAmplitudes)
    i->setModified(b);
}

std::shared_ptr<ComPWA::FunctionTree> SequentialPartialAmplitude::tree(
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

void SequentialPartialAmplitude::parameters(ParameterList &list) {
  Amplitude::parameters(list);
  for (auto i : PartialAmplitudes) {
    i->parameters(list);
  }
}

void SequentialPartialAmplitude::updateParameters(const ParameterList &list) {
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
