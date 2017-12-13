// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.
//

#include <numeric>
#include "Physics/IncoherentIntensity.hpp"

using namespace ComPWA::Physics;

std::shared_ptr<IncoherentIntensity>
IncoherentIntensity::Factory(std::shared_ptr<PartList> partL,
                             std::shared_ptr<Kinematics> kin,
                             const boost::property_tree::ptree &pt) {
  LOG(trace) << " IncoherentIntensity::Factory() | Construction....";

  auto obj = std::make_shared<IncoherentIntensity>();

  // Name is not required - default value 'empty'
  obj->Name = (pt.get<std::string>("<xmlattr>.Name", "empty"));

  std::shared_ptr<DoubleParameter> strength;
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter") {
      // Parameter (e.g. Mass)
      if (v.second.get<std::string>("<xmlattr>.Type") != "Strength")
        continue;
      auto tmp = DoubleParameter();
      tmp.load(v.second);
      strength = std::make_shared<DoubleParameter>(tmp);
    } else if (v.first == "CoherentIntensity") {
      obj->addIntensity(
          ComPWA::Physics::CoherentIntensity::Factory(partL, kin, v.second));
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
IncoherentIntensity::Save(std::shared_ptr<IncoherentIntensity> obj) {

  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->name());
  pt.add_child("Parameter", obj->Strength->save());
  pt.put("Parameter.<xmlattr>.Type", "Strength");
  for (auto i : obj->intensities()) {
    // TODO: we have to implement a memeber function Save() in AmpIntensity
    // interface and use it here
    auto ptr = std::dynamic_pointer_cast<CoherentIntensity>(i);
    pt.add_child("CoherentIntensity", CoherentIntensity::Save(ptr));
  }
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
  for (int i = 0; i < Intensities.size(); i++) {
    std::vector<double> params;
    Intensities.at(i)->parametersFast(params);
    if (parameters.at(i) != params) { // recalculate normalization
      parameters.at(i) = params;
      normValues.at(i) =
          1 / (Tools::Integral(Intensities.at(i), PhspSample, PhspVolume));
      normValues.at(i) *= Intensities.at(i)->strength();
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
    LOG(error) << "IncoherentIntensity::GetComponent() | You're requesting the "
                  "full object! So just copy it!";
    return std::shared_ptr<AmpIntensity>();
  }

  bool found = false;
  // Do we want to have a combination of CoherentIntensities?
  std::vector<std::string> names = splitString(name);
  auto icIn = std::shared_ptr<AmpIntensity>(this->clone(name));
  icIn->reset();
  for (auto i : names) {
    for (int j = 0; j < Intensities.size(); j++) {
      if (i == Intensities.at(j)->name()) {
        std::dynamic_pointer_cast<IncoherentIntensity>(icIn)->addIntensity(
            Intensities.at(j));
        found = true;
      }
    }
  }
  // Did we find something?
  if (found)
    return icIn;

  // Search for components in subsequent intensities
  for (auto i : Intensities) {
    try {
      auto r = i->component(name);
      std::dynamic_pointer_cast<IncoherentIntensity>(icIn)->addIntensity(r);
      found = true;
    } catch (std::exception &ex) {
    }
  }

  // Nothing found
  if (!found) {
    throw std::runtime_error(
        "InCoherentIntensity::GetComponent() | Component " + name +
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

  size_t sampleSize = sample.mDoubleValue(0)->values().size();


  auto tr = std::make_shared<FunctionTree>(
      "IncoherentIntens(" + name() + ")" + suffix, MDouble("",sampleSize),
      std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  tr->createLeaf("Strength", Strength,
                 "IncoherentIntens(" + name() + ")" + suffix);
  tr->createNode("SumOfCoherentIntens", MDouble("",sampleSize),
                 std::make_shared<AddAll>(ParType::MDOUBLE),
                 "IncoherentIntens(" + name() + ")" + suffix);
  for (auto i : Intensities) {
    tr->insertTree(i->tree(kin, sample, phspSample, toySample, nEvtVar),
                   "SumOfCoherentIntens");
  }
  return tr;
}

void IncoherentIntensity::parameters(ComPWA::ParameterList &list) {
  list.addParameter(Strength);
  for (auto i : Intensities) {
    i->parameters(list);
  }
}

void IncoherentIntensity::updateParameters(const ParameterList &list) {
  std::shared_ptr<DoubleParameter> p;
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
