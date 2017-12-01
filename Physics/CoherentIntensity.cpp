// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/Efficiency.hpp"
#include "Physics/CoherentIntensity.hpp"

using namespace ComPWA::Physics;

double CoherentIntensity::intensity(const DataPoint &point) const {
  std::complex<double> result(0., 0.);
  for (auto i : Amplitudes)
    result += i->evaluate(point);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return strength() * std::norm(result) * Eff->evaluate(point);
};

std::shared_ptr<CoherentIntensity>
CoherentIntensity::Factory(std::shared_ptr<PartList> partL,
                           std::shared_ptr<Kinematics> kin,
                           const boost::property_tree::ptree &pt) {
  LOG(trace) << " CoherentIntensity::Factory() | Construction....";
  auto obj = std::make_shared<CoherentIntensity>();
  obj->Name = (pt.get<std::string>("<xmlattr>.Name"));

  //  boost::property_tree::xml_writer_settings<char> settings('\t', 1);
  //  write_xml(std::cout,pt);

  std::shared_ptr<DoubleParameter> strength;
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter") {
      // Parameter (e.g. Mass)
      if (v.second.get<std::string>("<xmlattr>.Type") != "Strength")
        continue;
      auto tmp = DoubleParameter();
      tmp.load(v.second);
      strength = std::make_shared<DoubleParameter>(tmp);
    } else if (v.first == "Amplitude" &&
               v.second.get<std::string>("<xmlattr>.Class") ==
                   "SequentialPartialAmplitude") {
      obj->addAmplitude(
          ComPWA::Physics::HelicityFormalism::SequentialPartialAmplitude::
              Factory(partL, kin, v.second));
    } else {
      // ignored further settings. Should we throw an error?
    }
  }

  if (strength)
    obj->Strength = strength;
  else {
    obj->Strength = (std::make_shared<ComPWA::DoubleParameter>("", 1.0));
    obj->Strength->fixParameter(true);
  }

  obj->setPhspVolume(kin->phspVolume());

  return obj;
}

boost::property_tree::ptree
CoherentIntensity::Save(std::shared_ptr<CoherentIntensity> obj) {

  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->name());
  pt.add_child("Parameter", obj->Strength->save());
  pt.put("Parameter.<xmlattr>.Type", "Strength");

  for (auto i : obj->amplitudes()) {
    pt.add_child("Amplitude",
                 HelicityFormalism::SequentialPartialAmplitude::Save(i));
  }
  return pt;
}

std::shared_ptr<ComPWA::AmpIntensity>
CoherentIntensity::component(std::string name) {

  // The whole object?
  if (name == Name) {
    LOG(error) << "CoherentIntensity::GetComponent() | You're requesting the "
                  "full object! So just copy it!";
    return std::shared_ptr<AmpIntensity>();
  }

  bool found = false;
  // Do we want to have a combination of coherentintensities?
  std::vector<std::string> splitNames = splitString(name);
  auto icIn = std::shared_ptr<AmpIntensity>(this->clone(name));
  icIn->reset(); // delete all existing amplitudes
  for (auto i : splitNames) {
    for (auto j : Amplitudes) {
      if (i == j->name()) {
        std::dynamic_pointer_cast<CoherentIntensity>(icIn)->addAmplitude(j);
        found = true;
      }
    }
  }

  // Nothing found
  if (!found) {
    throw std::runtime_error(
        "CoherentIntensity::GetComponent() | Component " + name +
        " could not be found in CoherentIntensity " + Name + ".");
  }

  return icIn;
}

//! Getter function for basic amp tree
std::shared_ptr<ComPWA::FunctionTree>
CoherentIntensity::tree(std::shared_ptr<Kinematics> kin,
                           const ComPWA::ParameterList &sample,
                           const ComPWA::ParameterList &phspSample,
                           const ComPWA::ParameterList &toySample,
                           unsigned int nEvtVar, std::string suffix) {

  unsigned int effId = nEvtVar;
  unsigned int weightId = nEvtVar + 1;
  int phspSampleSize = phspSample.GetMultiDouble(0)->numValues();

  std::shared_ptr<MultiDouble> weightPhsp = phspSample.GetMultiDouble(weightId);
  double sumWeights =
      std::accumulate(weightPhsp->first(), weightPhsp->last(), 0.0);
  std::shared_ptr<MultiDouble> eff = phspSample.GetMultiDouble(effId);

  std::shared_ptr<FunctionTree> tr(new FunctionTree());

  tr->createHead("CoherentIntensity(" + name() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  tr->createLeaf("Strength", Strength,
                 "CoherentIntensity(" + name() + ")" + suffix);
  tr->createNode("SumSquared",
                 std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                 "CoherentIntensity(" + name() + ")" + suffix);
  tr->insertTree(setupBasicTree(kin, sample, toySample), "SumSquared");

  // Normalization
  // create a new FunctionTree to make sure that nodes with the same name do
  // not interfere
  std::shared_ptr<FunctionTree> trNorm(new FunctionTree());
  trNorm->createHead("Normalization(" + name() + ")" + suffix,
                     std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)));
  trNorm->createNode("Integral",
                     std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
                     "Normalization(" + name() + ")" + suffix);
  trNorm->createLeaf("PhspVolume", kin->phspVolume(), "Integral");
  trNorm->createLeaf("InverseSampleWeights", 1 / ((double)sumWeights),
                     "Integral");
  trNorm->createNode("Sum",
                     std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                     "Integral");
  trNorm->createNode("IntensityWeighted",
                     std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)),
                     "Sum", phspSampleSize, false);
  trNorm->createLeaf("Efficiency", eff, "IntensityWeighted");
  trNorm->createLeaf("EventWeight", weightPhsp, "IntensityWeighted");
  trNorm->createNode("Intensity",
                     std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                     "IntensityWeighted", phspSampleSize,
                     false); //|T_{ev}|^2
  trNorm->insertTree(setupBasicTree(kin, phspSample, toySample, "_norm"),
                     "Intensity");

  tr->insertTree(trNorm, "CoherentIntensity(" + name() + ")" + suffix);
  return tr;
}

std::shared_ptr<ComPWA::FunctionTree> CoherentIntensity::setupBasicTree(
    std::shared_ptr<Kinematics> kin, const ParameterList &sample,
    const ParameterList &phspSample, std::string suffix) const {

  int sampleSize = sample.GetMultiDouble(0)->numValues();
  int phspSampleSize = phspSample.GetMultiDouble(0)->numValues();

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

  newTree->createHead("SumOfAmplitudes" + suffix, maddStrat, sampleSize);

  for (auto i : Amplitudes) {
    std::shared_ptr<ComPWA::FunctionTree> resTree =
        i->tree(kin, sample, phspSample, "");
    if (!resTree->sanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Resonance tree didn't pass sanity check!");
    resTree->recalculate();
    newTree->insertTree(resTree, "SumOfAmplitudes" + suffix);
  }

  LOG(debug) << "AmpSumIntensity::setupBasicTree(): tree constructed!!";
  return newTree;
}

void CoherentIntensity::parameters(ComPWA::ParameterList &list) {
  list.AddParameter(Strength);
  for (auto i : Amplitudes) {
    i->parameters(list);
  }
}

void CoherentIntensity::updateParameters(const ParameterList &list) {
  std::shared_ptr<DoubleParameter> p;
  try {
    p = list.GetDoubleParameter(Strength->name());
  } catch (std::exception &ex) {
  }
  if (p)
    Strength->updateParameter(p);

  for (auto i : Amplitudes)
    i->updateParameters(list);

  return;
}
