// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/Efficiency.hpp"
#include "Physics/CoherentIntensity.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

double CoherentIntensity::Intensity(const dataPoint &point) const {
  std::complex<double> result(0., 0.);
  for (auto i : _seqDecays)
    result += i->Evaluate(point);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return Strength() * std::norm(result) * _eff->Evaluate(point);
};

std::shared_ptr<CoherentIntensity>
CoherentIntensity::Factory(std::shared_ptr<PartList> partL,
                           std::shared_ptr<Kinematics> kin,
                           const boost::property_tree::ptree &pt) {
  LOG(trace) << " CoherentIntensity::Factory() | Construction....";
  auto obj = std::make_shared<CoherentIntensity>();
  obj->_name = (pt.get<std::string>("<xmlattr>.Name"));

  //  boost::property_tree::xml_writer_settings<char> settings('\t', 1);
  //  write_xml(std::cout,pt);

  std::shared_ptr<DoubleParameter> strength;
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter") {
      // Parameter (e.g. Mass)
      if (v.second.get<std::string>("<xmlattr>.Type") != "Strength")
        continue;
      auto tmp = ComPWA::DoubleParameterFactory(v.second);
      strength = std::make_shared<DoubleParameter>(tmp);
    } else if (v.first == "Amplitude" &&
               v.second.get<std::string>("<xmlattr>.Class") ==
                   "SequentialPartialAmplitude") {
      obj->AddAmplitude(
          ComPWA::Physics::HelicityFormalism::SequentialPartialAmplitude::
              Factory(partL, kin, v.second));
    } else {
      // ignored further settings. Should we throw an error?
    }
  }

  if (strength)
    obj->_strength = strength;
  else {
    obj->_strength = (std::make_shared<ComPWA::DoubleParameter>("", 1.0));
    obj->_strength->SetParameterFixed();
  }

  obj->SetPhspVolume(kin->GetPhspVolume());

  return obj;
}

boost::property_tree::ptree
CoherentIntensity::Save(std::shared_ptr<CoherentIntensity> obj) {

  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->Name());
  pt.add_child("Parameter", ComPWA::DoubleParameterSave(*obj->_strength.get()));
  pt.put("Parameter.<xmlattr>.Type", "Strength");

  for (auto i : obj->GetAmplitudes()) {
    pt.add_child("Amplitude", SequentialPartialAmplitude::Save(i));
  }
  return pt;
}

std::shared_ptr<AmpIntensity>
CoherentIntensity::GetComponent(std::string name) {

  // The whole object?
  if (name == Name()) {
    LOG(error) << "CoherentIntensity::GetComponent() | You're requesting the "
                  "full object! So just copy it!";
    return std::shared_ptr<AmpIntensity>();
  }

  bool found = false;
  // Do we want to have a combination of coherentintensities?
  std::vector<std::string> splitNames = splitString(name);
  auto icIn = std::shared_ptr<AmpIntensity>(this->Clone(name));
  icIn->Reset(); // delete all existing amplitudes
  for (auto i : splitNames) {
    for (auto j : _seqDecays) {
      if (i == j->GetName()) {
        std::dynamic_pointer_cast<CoherentIntensity>(icIn)->AddAmplitude(j);
        found = true;
      }
    }
  }

  // Nothing found
  if (!found) {
    throw std::runtime_error(
        "CoherentIntensity::GetComponent() | Component " + name +
        " could not be found in CoherentIntensity " + Name() + ".");
  }

  return icIn;
}

//! Getter function for basic amp tree
std::shared_ptr<ComPWA::FunctionTree>
CoherentIntensity::GetTree(std::shared_ptr<Kinematics> kin,
                           const ComPWA::ParameterList &sample,
                           const ComPWA::ParameterList &phspSample,
                           const ComPWA::ParameterList &toySample,
                           unsigned int nEvtVar, std::string suffix) {

  unsigned int effId = nEvtVar;
  unsigned int weightId = nEvtVar + 1;
  int phspSampleSize = phspSample.GetMultiDouble(0)->GetNValues();

  std::shared_ptr<MultiDouble> weightPhsp = phspSample.GetMultiDouble(weightId);
  double sumWeights =
      std::accumulate(weightPhsp->Begin(), weightPhsp->End(), 0.0);
  std::shared_ptr<MultiDouble> eff = phspSample.GetMultiDouble(effId);

  std::shared_ptr<FunctionTree> tr(new FunctionTree());

  tr->CreateHead("CoherentIntensity(" + Name() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  tr->CreateLeaf("Strength", _strength,
                 "CoherentIntensity(" + Name() + ")" + suffix);
  tr->CreateNode("SumSquared",
                 std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                 "CoherentIntensity(" + Name() + ")" + suffix);
  tr->InsertTree(setupBasicTree(kin, sample, toySample), "SumSquared");

  // Normalization
  // create a new FunctionTree to make sure that nodes with the same name do
  // not interfere
  std::shared_ptr<FunctionTree> trNorm(new FunctionTree());
  trNorm->CreateHead("Normalization(" + Name() + ")" + suffix,
                     std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)));
  trNorm->CreateNode("Integral",
                     std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
                     "Normalization(" + Name() + ")" + suffix);
  trNorm->CreateLeaf("PhspVolume", kin->GetPhspVolume(), "Integral");
  trNorm->CreateLeaf("InverseSampleWeights", 1 / ((double)sumWeights),
                     "Integral");
  trNorm->CreateNode("Sum",
                     std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                     "Integral");
  trNorm->CreateNode("IntensityWeighted",
                     std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)),
                     "Sum", phspSampleSize, false);
  trNorm->CreateLeaf("Efficiency", eff, "IntensityWeighted");
  trNorm->CreateLeaf("EventWeight", weightPhsp, "IntensityWeighted");
  trNorm->CreateNode("Intensity",
                     std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                     "IntensityWeighted", phspSampleSize,
                     false); //|T_{ev}|^2
  trNorm->InsertTree(setupBasicTree(kin, phspSample, toySample, "_norm"),
                     "Intensity");

  tr->InsertTree(trNorm, "CoherentIntensity(" + Name() + ")" + suffix);
  return tr;
}

std::shared_ptr<FunctionTree> CoherentIntensity::setupBasicTree(
    std::shared_ptr<Kinematics> kin, const ParameterList &sample,
    const ParameterList &phspSample, std::string suffix) const {

  int sampleSize = sample.GetMultiDouble(0)->GetNValues();
  int phspSampleSize = phspSample.GetMultiDouble(0)->GetNValues();

  if (sampleSize == 0) {
    LOG(error) << "AmpSumIntensity::setupBasicTree() | "
                  "Data sample empty!";
    return std::shared_ptr<FunctionTree>();
  }
  if (phspSampleSize == 0) {
    LOG(error) << "AmpSumIntensity::setupBasicTree() | "
                  "Phsp sample empty!";
    return std::shared_ptr<FunctionTree>();
  }

  //------------ Setup FunctionTree ---------------------
  std::shared_ptr<FunctionTree> newTree(new FunctionTree());

  // Strategies needed
  std::shared_ptr<AddAll> maddStrat(new AddAll(ParType::MCOMPLEX));

  newTree->CreateHead("SumOfAmplitudes" + suffix, maddStrat, sampleSize);

  for (auto i : _seqDecays) {
    std::shared_ptr<FunctionTree> resTree =
        i->GetTree(kin, sample, phspSample, "");
    if (!resTree->SanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Resonance tree didn't pass sanity check!");
    resTree->Recalculate();
    newTree->InsertTree(resTree, "SumOfAmplitudes" + suffix);
  }

  LOG(debug) << "AmpSumIntensity::setupBasicTree(): tree constructed!!";
  return newTree;
}

void CoherentIntensity::GetParameters(ComPWA::ParameterList &list) {
  list.AddParameter(_strength);
  for (auto i : _seqDecays) {
    i->GetParameters(list);
  }
}

void CoherentIntensity::UpdateParameters(const ParameterList &list) {
  std::shared_ptr<DoubleParameter> p;
  try {
    p = list.GetDoubleParameter(_strength->name());
  } catch (std::exception &ex) {
  }
  if (p)
    _strength->UpdateParameter(p);

  for (auto i : _seqDecays)
    i->UpdateParameters(list);

  return;
}

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA
