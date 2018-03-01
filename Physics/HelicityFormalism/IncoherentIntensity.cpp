// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.
//

#include <numeric>
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::shared_ptr<IncoherentIntensity>
IncoherentIntensity::Factory(std::shared_ptr<PartList> partL,
                             std::shared_ptr<Kinematics> kin,
                             const boost::property_tree::ptree &pt) {
  LOG(trace) << " IncoherentIntensity::Factory() | Construction....";

  auto obj = std::make_shared<IncoherentIntensity>();

  // Name is not required - default value 'empty'
  obj->_name = (pt.get<std::string>("<xmlattr>.Name", "empty"));

  std::shared_ptr<DoubleParameter> strength;
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter") {
      // Parameter (e.g. Mass)
      if (v.second.get<std::string>("<xmlattr>.Type") != "Strength")
        continue;
      auto tmp = ComPWA::DoubleParameterFactory(v.second);
      strength = std::make_shared<DoubleParameter>(tmp);
    } else if (v.first == "CoherentIntensity") {
      obj->AddIntensity(
          ComPWA::Physics::HelicityFormalism::CoherentIntensity::Factory(
              partL, kin, v.second));
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
IncoherentIntensity::Save(std::shared_ptr<IncoherentIntensity> obj) {

  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->Name());
  pt.add_child("Parameter", ComPWA::DoubleParameterSave(*obj->_strength.get()));
  pt.put("Parameter.<xmlattr>.Type", "Strength");
  for (auto i : obj->GetIntensities()) {
    // TODO: we have to implement a memeber function Save() in AmpIntensity
    // interface and use it here
    auto ptr = std::dynamic_pointer_cast<CoherentIntensity>(i);
    pt.add_child("CoherentIntensity", CoherentIntensity::Save(ptr));
  }
  return pt;
}

double IncoherentIntensity::Intensity(const ComPWA::dataPoint &point) const {

  // We have to get around the constness of the interface definition.
  std::vector<std::vector<double>> parameters(_parameters);

  std::vector<double> normValues(_normValues);

  if (_intens.size() != parameters.size())
    parameters = std::vector<std::vector<double>>(_intens.size());

  if (_intens.size() != normValues.size())
    normValues = std::vector<double>(_intens.size());

  double result = 0;
  for (int i = 0; i < _intens.size(); i++) {
    std::vector<double> params;
    _intens.at(i)->GetParametersFast(params);
    if (parameters.at(i) != params) { // recalculate normalization
      parameters.at(i) = params;
      normValues.at(i) =
          1 / (Tools::Integral(_intens.at(i), _phspSample, phspVolume_));
      normValues.at(i) *= _intens.at(i)->Strength();
    }
    result += _intens.at(i)->Intensity(point) * normValues.at(i);
  }

  const_cast<std::vector<std::vector<double>> &>(_parameters) = parameters;
  const_cast<std::vector<double> &>(_normValues) = normValues;

  assert(!std::isnan(result) &&
         "IncoherentIntensity::Intensity() | Result is NaN!");
  assert(!std::isinf(result) &&
         "IncoherentIntensity::Intensity() | Result is inf!");

  return (Strength() * result);
}

std::shared_ptr<AmpIntensity>
IncoherentIntensity::GetComponent(std::string name) {

  // The whole object?
  if (name == Name()) {
    LOG(error) << "IncoherentIntensity::GetComponent() | You're requesting the "
                  "full object! So just copy it!";
    return std::shared_ptr<AmpIntensity>();
  }

  bool found = false;
  // Do we want to have a combination of CoherentIntensities?
  std::vector<std::string> names = splitString(name);
  auto icIn = std::shared_ptr<AmpIntensity>(this->Clone(name));
  icIn->Reset();
  for (auto i : names) {
    for (int j = 0; j < _intens.size(); j++) {
      if (name == _intens.at(j)->Name()) {
        std::dynamic_pointer_cast<IncoherentIntensity>(icIn)->AddIntensity(
            _intens.at(j));
        if (names.size() == 1)
          return _intens.at(j);
        found = true;
      }
    }
  }

  // Did we find something?
  if (found)
    return icIn;

  // Search for components in subsequent intensities
  for (auto i : _intens) {
    try {
      icIn = i->GetComponent(name);
      found = true;
    } catch (std::exception &ex) {
    }
  }

  // Nothing found
  if (!found) {
    throw std::runtime_error(
        "InCoherentIntensity::GetComponent() | Component " + name +
        " could not be found in IncoherentIntensity " + Name() + ".");
  }

  return icIn;
}

std::shared_ptr<ComPWA::FunctionTree>
IncoherentIntensity::GetTree(std::shared_ptr<Kinematics> kin,
                             const ComPWA::ParameterList &sample,
                             const ComPWA::ParameterList &phspSample,
                             const ComPWA::ParameterList &toySample,
                             unsigned int nEvtVar, std::string suffix) {

  std::shared_ptr<FunctionTree> tr(new FunctionTree());

  tr->CreateHead("IncoherentIntens(" + Name() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  tr->CreateLeaf("Strength", _strength,
                 "IncoherentIntens(" + Name() + ")" + suffix);
  tr->CreateNode("SumOfCoherentIntens",
                 std::shared_ptr<Strategy>(new AddAll(ParType::MDOUBLE)),
                 "IncoherentIntens(" + Name() + ")" + suffix);
  for (auto i : _intens) {
    tr->InsertTree(i->GetTree(kin, sample, phspSample, toySample, nEvtVar),
                   "SumOfCoherentIntens");
  }
  return tr;
}

void IncoherentIntensity::GetParameters(ComPWA::ParameterList &list) {
  std::shared_ptr<DoubleParameter> tmp;
  try { // catch BadParameter
    tmp = list.GetDoubleParameter(_strength->GetName());
    // catch and throw std::runtime_error due to failed parameter comparisson
    try {
      if (*tmp == *_strength)
        _strength = tmp;
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.AddParameter(_strength);
  }

  for (auto i : GetIntensities()) {
    i->GetParameters(list);
  }
}

void IncoherentIntensity::UpdateParameters(const ParameterList &list) {
  std::shared_ptr<DoubleParameter> p;
  try {
    p = list.GetDoubleParameter(_strength->GetName());
  } catch (std::exception &ex) {
  }
  if (p)
    _strength->UpdateParameter(p);
  for (auto i : _intens)
    i->UpdateParameters(list);

  return;
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
