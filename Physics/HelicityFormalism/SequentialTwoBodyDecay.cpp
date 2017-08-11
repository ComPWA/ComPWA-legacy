// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <memory>
#include "Physics/HelicityFormalism/SequentialTwoBodyDecay.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::shared_ptr<ComPWA::Physics::Amplitude>
SequentialTwoBodyDecay::Factory(std::shared_ptr<Kinematics> kin,
                                const boost::property_tree::ptree &pt) {
  LOG(trace) << " SequentialTwoBodyDecay::Factory() | Construction....";
  auto obj = std::make_shared<SequentialTwoBodyDecay>();
  obj->SetName(pt.get<std::string>("<xmlattr>.Name", "empty"));

  auto mag = ComPWA::DoubleParameterFactory(pt.get_child("Magnitude"));
  obj->SetMagnitudeParameter(std::make_shared<DoubleParameter>(mag));
  auto phase = ComPWA::DoubleParameterFactory(pt.get_child("Phase"));
  obj->SetPhaseParameter(std::make_shared<DoubleParameter>(phase));

  auto prefactor = pt.get_child_optional("PreFactor");
  if (prefactor) {
    double r = prefactor.get().get<double>("<xmlattr>.Magnitude");
    double p = prefactor.get().get<double>("<xmlattr>.Phase");
    auto pref = std::polar(r, p);
    obj->SetPreFactor(pref);
  }

  for (const auto &v : pt.get_child("")) {
    if (v.first == "Resonance")
      obj->Add(PartialDecay::Factory(kin, v.second));
  }

  return std::static_pointer_cast<ComPWA::Physics::Amplitude>(obj);
}

boost::property_tree::ptree
SequentialTwoBodyDecay::Save(std::shared_ptr<ComPWA::Physics::Amplitude> amp) {

  auto obj = std::static_pointer_cast<SequentialTwoBodyDecay>(amp);
  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->GetName());
  pt.add_child("Magnitude", ComPWA::DoubleParameterSave(
                                *obj->GetMagnitudeParameter().get()));
  pt.add_child("Phase",
               ComPWA::DoubleParameterSave(*obj->GetPhaseParameter().get()));

  auto pref = amp->GetPreFactor();
  if (pref != std::complex<double>(1, 0)) {
    boost::property_tree::ptree ppp;
    ppp.put("<xmlattr>.Magnitude", std::abs(pref));
    ppp.put("<xmlattr>.Phase", std::arg(pref));
    pt.add_child("PreFactor", ppp);
  }

  for (auto i : obj->GetDecays()) {
    pt.add_child("Resonance", PartialDecay::Save(i));
  }
  return pt;
}

/**! Setup function tree */
std::shared_ptr<FunctionTree> SequentialTwoBodyDecay::GetTree(
    std::shared_ptr<Kinematics> kin, const ParameterList &sample,
    const ParameterList &toySample, std::string suffix) {

  std::shared_ptr<FunctionTree> tr(new FunctionTree());
  tr->createHead("Amplitude(" + GetName() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)));
  tr->createNode("Strength",
                 std::shared_ptr<Strategy>(new Complexify(ParType::COMPLEX)),
                 "Amplitude(" + GetName() + ")" + suffix);
  tr->createLeaf("Magnitude", _magnitude, "Strength");
  tr->createLeaf("Phase", _phase, "Strength");
  tr->createLeaf("PreFactor", GetPreFactor(),
                 "Amplitude(" + GetName() + ")" + suffix);

  for (auto i : _partDecays) {
    std::shared_ptr<FunctionTree> resTree =
        i->GetTree(kin, sample, toySample, "");
    if (!resTree->sanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Amplitude tree didn't pass sanity check!");
    resTree->recalculate();
    tr->insertTree(resTree, "Amplitude(" + GetName() + ")" + suffix);
  }

  return tr;
};

void SequentialTwoBodyDecay::GetParameters(ParameterList &list) {
  Amplitude::GetParameters(list);
  for (auto i : _partDecays) {
    i->GetParameters(list);
  }
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
