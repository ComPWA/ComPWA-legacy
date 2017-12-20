// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/Efficiency.hpp"
#include "Physics/CoherentIntensity.hpp"

using namespace ComPWA;
using namespace ComPWA::Physics;
using namespace ComPWA::Physics::HelicityFormalism;

CoherentIntensity::CoherentIntensity(std::shared_ptr<ComPWA::PartList> partL,
                  std::shared_ptr<ComPWA::Kinematics> kin,
                  const boost::property_tree::ptree &pt) : PhspVolume(1.0) {
  load(partL, kin, pt);
}

double CoherentIntensity::intensity(const DataPoint &point) const {
  std::complex<double> result(0., 0.);
  for (auto i : Amplitudes)
    result += i->evaluate(point);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return strength() * std::norm(result) * Eff->evaluate(point);
};

void CoherentIntensity::load(std::shared_ptr<PartList> partL,
                             std::shared_ptr<Kinematics> kin,
                             const boost::property_tree::ptree &pt) {
  LOG(trace) << "CoherentIntensity::load() | Construction....";
  if (pt.get<std::string>("<xmlattr>.Class") != "Coherent")
    throw BadConfig("CoherentIntensity::Factory() | Property tree seems to "
                    "not containt a configuration for an "
                    "CoherentIntensity!");
  Name = (pt.get<std::string>("<xmlattr>.Name"));
  Strength = std::make_shared<ComPWA::FitParameter>("Strength_"+Name, 1.0);
  setPhspVolume(kin->phspVolume());
  
  for (const auto &v : pt.get_child("")) {
    // Mass parameter
    if (v.first == "Parameter" &&
        v.second.get<std::string>("<xmlattr>.Type") != "Strength") {
      Strength = std::make_shared<FitParameter>(v.second);
    } else if (v.first == "Amplitude" &&
               v.second.get<std::string>("<xmlattr>.Class") ==
                   "SequentialPartialAmplitude") {
      addAmplitude(
          std::make_shared<SequentialPartialAmplitude>(partL, kin, v.second));
    } else {
      // ignored further settings. Should we throw an error?
    }
  }
}

boost::property_tree::ptree CoherentIntensity::save() const {
  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", name());
  pt.add_child("Parameter", Strength->save());
  pt.put("Parameter.<xmlattr>.Type", "Strength");

  for (auto i : Amplitudes)
    pt.add_child("Amplitude", i->save());
  
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
  auto icIn = std::make_shared<CoherentIntensity>(*this);
  icIn->setName(name);
  icIn->reset(); // delete all existing amplitudes
  for (auto i : splitNames) {
    for (auto j : Amplitudes) {
      if (i == j->name()) {
        icIn->addAmplitude(j);
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

  size_t n = sample.mDoubleValue(0)->values().size();
  size_t phspSize = phspSample.mDoubleValue(0)->values().size();

  // Efficiency values are stored on the next to last element of the
  // ParameterList
  std::shared_ptr<Value<std::vector<double>>> eff =
      phspSample.mDoubleValues().end()[-2];
  // Weights are stored on the last element of the ParameterList
  std::shared_ptr<Value<std::vector<double>>> weightPhsp =
      phspSample.mDoubleValues().end()[-1];

  double sumWeights = std::accumulate(weightPhsp->values().begin(),
                                      weightPhsp->values().end(), 0.0);

  auto tr = std::make_shared<FunctionTree>(
      "CoherentIntensity(" + name() + ")" + suffix, MDouble("", n),
      std::make_shared<MultAll>(ParType::MDOUBLE));

  tr->createLeaf("Strength", Strength,
                 "CoherentIntensity(" + name() + ")" + suffix);
  tr->createNode("SumSquared",
                 std::make_shared<AbsSquare>(ParType::MDOUBLE),
                 "CoherentIntensity(" + name() + ")" + suffix);
  tr->insertTree(setupBasicTree(kin, sample, toySample), "SumSquared");

  // Normalization
  // create a new FunctionTree to make sure that nodes with the same name do
  // not interfere
  auto trNorm = std::make_shared<FunctionTree>(
      "Normalization(" + name() + ")" + suffix,
      std::make_shared<Value<double>>(),
      std::make_shared<Inverse>(ParType::DOUBLE));
  trNorm->createNode("Integral",
                     std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
                     "Normalization(" + name() + ")" + suffix);
  trNorm->createLeaf("PhspVolume", kin->phspVolume(), "Integral");
  trNorm->createLeaf("InverseSampleWeights", 1 / ((double)sumWeights),
                     "Integral");
  trNorm->createNode("Sum",
                     std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                     "Integral");
  trNorm->createNode("IntensityWeighted", MDouble("", phspSize),
                     std::make_shared<MultAll>(ParType::MDOUBLE), "Sum");
  trNorm->createLeaf("Efficiency", eff, "IntensityWeighted");
  trNorm->createLeaf("EventWeight", weightPhsp, "IntensityWeighted");
  trNorm->createNode("Intensity", MDouble("", phspSize),
                     std::make_shared<AbsSquare>(ParType::MDOUBLE),
                     "IntensityWeighted"); //|T_{ev}|^2
  trNorm->insertTree(setupBasicTree(kin, phspSample, toySample, "_norm"),
                     "Intensity");

  tr->insertTree(trNorm, "CoherentIntensity(" + name() + ")" + suffix);
  return tr;
}

std::shared_ptr<ComPWA::FunctionTree> CoherentIntensity::setupBasicTree(
    std::shared_ptr<Kinematics> kin, const ParameterList &sample,
    const ParameterList &phspSample, std::string suffix) const {

  size_t n = sample.mDoubleValue(0)->values().size();

  if (n == 0) {
    LOG(error) << "AmpSumIntensity::setupBasicTree() | "
                  "Data sample empty!";
    return std::shared_ptr<FunctionTree>();
  }

  //------------ Setup FunctionTree ---------------------
  auto newTree = std::make_shared<FunctionTree>(
      "SumOfAmplitudes" + suffix, MComplex("", n),
      std::make_shared<AddAll>(ParType::MCOMPLEX));

  for (auto i : Amplitudes) {
    std::shared_ptr<ComPWA::FunctionTree> resTree =
        i->tree(kin, sample, phspSample, "");
    if (!resTree->sanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Resonance tree didn't pass sanity check!");
    resTree->parameter();
    newTree->insertTree(resTree, "SumOfAmplitudes" + suffix);
  }

  LOG(debug) << "AmpSumIntensity::setupBasicTree(): tree constructed!!";
  return newTree;
}

void CoherentIntensity::parameters(ComPWA::ParameterList &list) {
  Strength = list.addUniqueParameter(Strength);
  for (auto i : Amplitudes) {
    i->parameters(list);
  }
}

void CoherentIntensity::updateParameters(const ParameterList &list) {
  std::shared_ptr<FitParameter> p;
  try {
    p = FindParameter(Strength->name(), list);
  } catch (std::exception &ex) {
  }
  if (p)
    Strength->updateParameter(p);

  for (auto i : Amplitudes)
    i->updateParameters(list);

  return;
}
