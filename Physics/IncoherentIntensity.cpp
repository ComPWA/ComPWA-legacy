// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.
//

#include <numeric>
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/CoherentIntensity.hpp"

using namespace ComPWA::Physics;

IncoherentIntensity::IncoherentIntensity(
    std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
    const boost::property_tree::ptree &pt) {
  load(partL, kin, pt);
}

void IncoherentIntensity::load(std::shared_ptr<PartList> partL,
                               std::shared_ptr<Kinematics> kin,
                               const boost::property_tree::ptree &pt) {
  LOG(TRACE) << " IncoherentIntensity::Factory() | Construction....";
  if (pt.get<std::string>("<xmlattr>.Class") != "Incoherent")
    throw BadConfig("IncoherentIntensity::Factory() | Property tree seems to "
                    "not containt a configuration for an "
                    "IncoherentIntensity!");

  // Name is not required - default value 'empty'
  Name = (pt.get<std::string>("<xmlattr>.Name", "empty"));
  Strength = (std::make_shared<ComPWA::FitParameter>("Strength_" + Name, 1.0));
  setPhspVolume(kin->phspVolume());

  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter" &&
        v.second.get<std::string>("<xmlattr>.Type") == "Strength") {
      Strength = std::make_shared<FitParameter>(v.second);
    } else if (v.first == "Intensity" &&
               v.second.get<std::string>("<xmlattr>.Class") == "Coherent") {
      addIntensity(std::make_shared<CoherentIntensity>(partL, kin, v.second));
    } else {
      // ignored further settings. Should we throw an error?
    }
  }
}

boost::property_tree::ptree IncoherentIntensity::save() const {
  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", name());
  pt.add_child("Parameter", Strength->save());
  pt.put("Parameter.<xmlattr>.Type", "Strength");
  for (auto i : Intensities)
    pt.add_child("CoherentIntensity", i->save());

  return pt;
}

double IncoherentIntensity::intensity(const ComPWA::DataPoint &point) const {

  // We have to get around the constness of the interface definition.
  std::vector<std::vector<double>> parameters(Parameters);

  std::vector<double> normValues(NormalizationValues);

  if (Intensities.size() != parameters.size())
    parameters = std::vector<std::vector<double>>(Intensities.size());

  if (Intensities.size() != normValues.size())
    normValues = std::vector<double>(Intensities.size());

  double result = 0;
  for (unsigned int i = 0; i < Intensities.size(); ++i) {
    std::vector<double> params;
    Intensities.at(i)->parametersFast(params);
    if (parameters.at(i) != params) { // recalculate normalization
      parameters.at(i) = params;
      normValues.at(i) =
          1 / (Tools::Integral(Intensities.at(i), PhspSample, PhspVolume));
    }
    result += Intensities.at(i)->intensity(point) * normValues.at(i);
  }

  const_cast<std::vector<std::vector<double>> &>(Parameters) = parameters;
  const_cast<std::vector<double> &>(NormalizationValues) = normValues;

  assert(!std::isnan(result) &&
         "IncoherentIntensity::Intensity() | Result is NaN!");
  assert(!std::isinf(result) &&
         "IncoherentIntensity::Intensity() | Result is inf!");

  return (strength() * result);
}

std::shared_ptr<ComPWA::AmpIntensity>
IncoherentIntensity::component(std::string name) {

  // The whole object?
  if (name == Name) {
//    LOG(ERROR) << "IncoherentIntensity::GetComponent() | You're requesting the "
//                  "full object! So just copy it!";
//    return std::shared_ptr<AmpIntensity>();
    return shared_from_this();
  }
  
  // Do we want to have a combination of CoherentIntensities?
  std::vector<std::string> names = splitString(name);
  
  // In case the requested is a direct component we return a CoherentIntensity
  // object
  if (names.size() == 1) {
    for (unsigned int j = 0; j < Intensities.size(); ++j) {
      if (names.at(0) == Intensities.at(j)->name()) {
        return Intensities.at(j);
      }
    }
  }

  // Otherwise we hava multiple components and we build a IncoherentIntensity
  auto icIn = std::make_shared<IncoherentIntensity>(*this);
  icIn->setName(name);
  icIn->reset();

  for (unsigned int j = 0; j < Intensities.size(); ++j) {
    for (auto i : names) {
      if (i == Intensities.at(j)->name()) {
        icIn->addIntensity(Intensities.at(j));
      }
    }
  }

  // Nothing found
  if (names.size() != icIn->intensities().size()) {
    throw std::runtime_error(
        "InCoherentIntensity::component() | Component " + name +
        " could not be found in IncoherentIntensity " + Name + ".");
  }

  return icIn;
}

std::shared_ptr<ComPWA::AmpIntensity>
IncoherentIntensity::component(std::string name, std::string resName, 
    std::string daug1Name, std::string daug2Name, int L, int S) {
  std::string lStr = "L_" + std::to_string(L) + ".0";
  if (L < 0) lStr = "";
  std::string sStr = "S_" + std::to_string(S) + ".0";
  if (S < 0) sStr = "";
  std::string lsStr = lStr + "_" + sStr;
  if (name == "") {
    name = resName + "_to_" + daug1Name + "+" + daug2Name;
    if (L >= 0) name += "_" + lStr;
    if (S >= 0) name += "_" + sStr;
  }

  auto icIn = std::make_shared<IncoherentIntensity>(*this);
  icIn->setName(name);
  icIn->reset(); // delete all existing amplitudes

  bool found(false);
  for (auto ampIntensity : Intensities) {
    auto coherentIntensity = 
        std::dynamic_pointer_cast<ComPWA::Physics::CoherentIntensity>(ampIntensity);
    std::shared_ptr<ComPWA::AmpIntensity> comp;
    try {
      comp = coherentIntensity->component(name, resName, daug1Name, daug2Name, L, S);
    } catch (std::exception &ex) {
      continue;
    }
    found = true;
    icIn->addIntensity(comp);
  }

  // Nothing found
  if (!found) {
    throw std::runtime_error(
        "InCoherentIntensity::component() | Component " + name +
        " could not be found in IncoherentIntensity " + Name + ".");
  }

  return icIn;
}

std::shared_ptr<ComPWA::FunctionTree>
IncoherentIntensity::tree(std::shared_ptr<Kinematics> kin,
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
      "IncoherentIntens(" + name() + ")" + suffix, MDouble("", n),
      std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  tr->createLeaf("Strength", Strength,
                 "IncoherentIntens(" + name() + ")" + suffix);
  tr->createNode("SumOfCoherentIntens", MDouble("", n),
                 std::make_shared<AddAll>(ParType::MDOUBLE),
                 "IncoherentIntens(" + name() + ")" + suffix);

  for (auto i : Intensities) {
    std::string name = i->name();
    // Construct tree for normalization
    auto normTree = std::make_shared<FunctionTree>(
        "Normalization(" + name + ")" + suffix,
        std::make_shared<Value<double>>(),
        std::make_shared<Inverse>(ParType::DOUBLE));
    normTree->createNode(
        "Integral", std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
        "Normalization(" + name + ")" + suffix);
    normTree->createLeaf("PhspVolume", kin->phspVolume(), "Integral");
    normTree->createLeaf("InverseSampleWeights", 1 / ((double)sumWeights),
                         "Integral");
    normTree->createNode("Sum",
                         std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                         "Integral");
    normTree->createNode("IntensityWeighted", MDouble("", phspSize),
                         std::make_shared<MultAll>(ParType::MDOUBLE), "Sum");
    normTree->createLeaf("Efficiency", eff, "IntensityWeighted");
    normTree->createLeaf("EventWeight", weightPhsp, "IntensityWeighted");
    normTree->insertTree(
        i->tree(kin, phspSample, phspSample, toySample, nEvtVar),
        "IntensityWeighted");

    // Construct tree and add normalization tree
    auto intensTree = i->tree(kin, sample, phspSample, toySample, nEvtVar);
    intensTree->insertTree(normTree,
                           "CoherentIntensity(" + name + ")" + suffix);

    tr->insertTree(intensTree, "SumOfCoherentIntens");
  }
  return tr;
}

void IncoherentIntensity::parameters(ComPWA::ParameterList &list) {
  Strength = list.addUniqueParameter(Strength);
  for (auto i : Intensities) {
    i->parameters(list);
  }
}

void IncoherentIntensity::updateParameters(const ParameterList &list) {
  std::shared_ptr<FitParameter> p;
  try {
    p = FindParameter(Strength->name(), list);
  } catch (std::exception &ex) {
  }
  if (p)
    Strength->updateParameter(p);
  for (auto i : Intensities)
    i->updateParameters(list);

  return;
}
