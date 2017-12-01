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
  obj->SetName(pt.get<std::string>("<xmlattr>.Name", "empty"));

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
    } else if (v.first == "PartialAmplitude" && v.second.get<std::string>("<xmlattr>.Class") == "HelicityDecay") {
      obj->Add(HelicityDecay::Factory(partL, kin, v.second));
    } else if (v.first == "PreFactor") {
      double r = v.second.get<double>("<xmlattr>.Magnitude");
      double p = v.second.get<double>("<xmlattr>.Phase");
      pref = std::polar(r, p);
    } else {
      // ignored further settings. Should we throw an error?
    }
  }


  if (mag)
    obj->SetMagnitudeParameter(mag);
  else
    throw BadParameter(
        "SequentialPartialAmplitude::Factory() | No magnitude parameter found.");

  if (phase)
    obj->SetPhaseParameter(phase);
  else
    throw BadParameter(
        "SequentialPartialAmplitude::Factory() | No phase parameter found.");

  obj->SetPreFactor(pref);

  return std::static_pointer_cast<ComPWA::Physics::Amplitude>(obj);
}

boost::property_tree::ptree
SequentialPartialAmplitude::Save(std::shared_ptr<ComPWA::Physics::Amplitude> amp) {

  auto obj = std::static_pointer_cast<SequentialPartialAmplitude>(amp);
  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->GetName());

  boost::property_tree::ptree tmp = obj->GetMagnitudeParameter()->save();
  tmp.put("<xmlattr>.Type", "Magnitude");
  pt.add_child("Parameter", tmp);

  tmp = obj->GetPhaseParameter()->save();
  tmp.put("<xmlattr>.Type", "Phase");
  pt.add_child("Parameter", tmp);

  auto pref = amp->GetPreFactor();
  if (pref != std::complex<double>(1, 0)) {
    boost::property_tree::ptree ppp;
    ppp.put("<xmlattr>.Magnitude", std::abs(pref));
    ppp.put("<xmlattr>.Phase", std::arg(pref));
    pt.add_child("PreFactor", ppp);
  }

  for (auto i : obj->GetDecays()) {
    pt.add_child("PartialAmplitude", HelicityDecay::Save(i));
  }
  return pt;
}

std::shared_ptr<ComPWA::FunctionTree> SequentialPartialAmplitude::GetTree(
    std::shared_ptr<Kinematics> kin, const ParameterList &sample,
    const ParameterList &toySample, std::string suffix) {

  std::shared_ptr<FunctionTree> tr(new FunctionTree());
  tr->CreateHead("Amplitude(" + GetName() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)));
  tr->CreateNode("Strength",
                 std::shared_ptr<Strategy>(new Complexify(ParType::COMPLEX)),
                 "Amplitude(" + GetName() + ")" + suffix);
  tr->CreateLeaf("Magnitude", _magnitude, "Strength");
  tr->CreateLeaf("Phase", _phase, "Strength");
  tr->CreateLeaf("PreFactor", GetPreFactor(),
                 "Amplitude(" + GetName() + ")" + suffix);

  for (auto i : _partDecays) {
    std::shared_ptr<FunctionTree> resTree =
        i->GetTree(kin, sample, toySample, "");
    if (!resTree->SanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Amplitude tree didn't pass sanity check!");
    resTree->Recalculate();
    tr->InsertTree(resTree, "Amplitude(" + GetName() + ")" + suffix);
  }

  return tr;
}

void SequentialPartialAmplitude::GetParameters(ParameterList &list) {
  Amplitude::GetParameters(list);
  for (auto i : _partDecays) {
    i->GetParameters(list);
  }
}

void SequentialPartialAmplitude::UpdateParameters(const ParameterList &par) {
  // Try to update magnitude
  std::shared_ptr<DoubleParameter> mag;
  try {
    mag = par.GetDoubleParameter(_magnitude->name());
  } catch (std::exception &ex) {
  }
  if (mag)
    _magnitude->updateParameter(mag);
  std::shared_ptr<DoubleParameter> phase;
  
  // Try to update phase
  try {
    phase = par.GetDoubleParameter(_phase->name());
  } catch (std::exception &ex) {
  }
  if (phase)
    _phase->updateParameter(phase);
  
  for (auto i : _partDecays)
    i->UpdateParameters(par);

  return;
}
