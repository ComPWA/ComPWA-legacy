// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <memory>
#include "Physics/SequentialPartialAmplitude.hpp"

using namespace ComPWA::Physics::HelicityFormalism;

std::shared_ptr<ComPWA::Physics::Amplitude>
SequentialPartialAmplitude::Factory(std::shared_ptr<PartList> partL,
                                    std::shared_ptr<Kinematics> kin,
                                    const boost::property_tree::ptree &pt) {
  LOG(trace) << " SequentialPartialAmplitude::Factory() | Construction....";
  auto obj = std::make_shared<SequentialPartialAmplitude>();
  obj->setName(pt.get<std::string>("<xmlattr>.Name", "empty"));

  std::shared_ptr<DoubleParameter> mag, phase;
  std::complex<double> pref(1, 0);
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter") {
      if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude") {
        auto tmp = DoubleParameter();
        tmp.load(v.second);
        mag = std::make_shared<DoubleParameter>(tmp);
      }
      if (v.second.get<std::string>("<xmlattr>.Type") == "Phase") {
        auto tmp = DoubleParameter();
        tmp.load(v.second);
        phase = std::make_shared<DoubleParameter>(tmp);
      }
    } else if (v.first == "PartialAmplitude" &&
               v.second.get<std::string>("<xmlattr>.Class") ==
                   "HelicityDecay") {
      obj->addPartialAmplitude(HelicityDecay::Factory(partL, kin, v.second));
    } else if (v.first == "PreFactor") {
      double r = v.second.get<double>("<xmlattr>.Magnitude");
      double p = v.second.get<double>("<xmlattr>.Phase");
      pref = std::polar(r, p);
    } else {
      // ignored further settings. Should we throw an error?
    }
  }

  if (mag)
    obj->setMagnitudeParameter(mag);
  else
    throw BadParameter("SequentialPartialAmplitude::Factory() | No magnitude "
                       "parameter found.");

  if (phase)
    obj->setPhaseParameter(phase);
  else
    throw BadParameter(
        "SequentialPartialAmplitude::Factory() | No phase parameter found.");

  obj->setPreFactor(pref);

  return std::static_pointer_cast<ComPWA::Physics::Amplitude>(obj);
}

boost::property_tree::ptree SequentialPartialAmplitude::Save(
    std::shared_ptr<ComPWA::Physics::Amplitude> amp) {

  auto obj = std::static_pointer_cast<SequentialPartialAmplitude>(amp);
  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->name());

  boost::property_tree::ptree tmp = obj->magnitudeParameter()->save();
  tmp.put("<xmlattr>.Type", "Magnitude");
  pt.add_child("Parameter", tmp);

  tmp = obj->phaseParameter()->save();
  tmp.put("<xmlattr>.Type", "Phase");
  pt.add_child("Parameter", tmp);

  auto pref = amp->preFactor();
  if (pref != std::complex<double>(1, 0)) {
    boost::property_tree::ptree ppp;
    ppp.put("<xmlattr>.Magnitude", std::abs(pref));
    ppp.put("<xmlattr>.Phase", std::arg(pref));
    pt.add_child("PreFactor", ppp);
  }

  for (auto i : obj->partialAmplitudes()) {
    pt.add_child("PartialAmplitude", HelicityDecay::Save(i));
  }
  return pt;
}

std::shared_ptr<ComPWA::FunctionTree> SequentialPartialAmplitude::tree(
    std::shared_ptr<Kinematics> kin, const ParameterList &sample,
    const ParameterList &toySample, std::string suffix) {

  std::shared_ptr<FunctionTree> tr(new FunctionTree());
  tr->createHead("Amplitude(" + name() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)));
  tr->createNode("Strength",
                 std::shared_ptr<Strategy>(new Complexify(ParType::COMPLEX)),
                 "Amplitude(" + name() + ")" + suffix);
  tr->createLeaf("Magnitude", Magnitude, "Strength");
  tr->createLeaf("Phase", Phase, "Strength");
  tr->createLeaf("PreFactor", preFactor(),
                 "Amplitude(" + name() + ")" + suffix);

  for (auto i : PartialAmplitudes) {
    std::shared_ptr<FunctionTree> resTree =
        i->tree(kin, sample, toySample, "");
    if (!resTree->sanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Amplitude tree didn't pass sanity check!");
    resTree->recalculate();
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

void SequentialPartialAmplitude::updateParameters(const ParameterList &par) {
  // Try to update magnitude
  std::shared_ptr<DoubleParameter> mag;
  try {
    mag = par.GetDoubleParameter(Magnitude->name());
  } catch (std::exception &ex) {
  }
  if (mag)
    Magnitude->updateParameter(mag);
  std::shared_ptr<DoubleParameter> phase;

  // Try to update phase
  try {
    phase = par.GetDoubleParameter(Phase->name());
  } catch (std::exception &ex) {
  }
  if (phase)
    Phase->updateParameter(phase);

  for (auto i : PartialAmplitudes)
    i->updateParameters(par);

  return;
}
