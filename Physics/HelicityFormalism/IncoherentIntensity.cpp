
//
//  IncoherentIntensity.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 21/02/2017.
//
//

#include <numeric>
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::shared_ptr<IncoherentIntensity>
IncoherentIntensity::Factory(const boost::property_tree::ptree &pt) {
  LOG(trace) << " IncoherentIntensity::Factory() | Construction....";

  auto obj = std::make_shared<IncoherentIntensity>();

  // Name is not required - default value 'empty'
  obj->_name=(pt.get<std::string>("<xmlattr>.Name", "empty"));

  auto ptCh = pt.get_child_optional("Strength");
  if (ptCh) {
    auto strength = ComPWA::DoubleParameterFactory(ptCh.get());
    obj->_strength=(std::make_shared<DoubleParameter>(strength));
  } else {
    obj->_strength=
        (std::make_shared<ComPWA::DoubleParameter>("", 1.0));
    obj->_strength->SetParameterFixed();
  }

  for (const auto &v : pt.get_child("")) {
    if (v.first == "CoherentIntensity")
      obj->AddIntensity(
          ComPWA::Physics::HelicityFormalism::CoherentIntensity::Factory(
              v.second));
  }
  return obj;
}

boost::property_tree::ptree
IncoherentIntensity::Save(std::shared_ptr<IncoherentIntensity> obj) {

  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->Name());
  pt.add_child("Strength",
               ComPWA::DoubleParameterSave(*obj->_strength.get()));
  for (auto i : obj->GetIntensities()) {
    // TODO: we have to implement a memeber function Save() in AmpIntensity
    // interface and use it here
    auto ptr = std::dynamic_pointer_cast<CoherentIntensity>(i);
    pt.add_child("CoherentIntensity", CoherentIntensity::Save(ptr));
  }
  return pt;
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
        if( names.size() == 1 ) return _intens.at(j);
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
IncoherentIntensity::GetTree(const ComPWA::ParameterList &sample,
                             const ComPWA::ParameterList &phspSample,
                             const ComPWA::ParameterList &toySample,
                             std::string suffix) {

  std::shared_ptr<FunctionTree> tr(new FunctionTree());

  tr->createHead("IncoherentIntens(" + Name() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  tr->createLeaf("Strength", _strength,
                 "IncoherentIntens(" + Name() + ")" + suffix);
  tr->createNode("SumOfCoherentIntens",
                 std::shared_ptr<Strategy>(new AddAll(ParType::MDOUBLE)),
                 "IncoherentIntens(" + Name() + ")" + suffix);
  for (auto i : _intens) {
    tr->insertTree(i->GetTree(sample, phspSample, toySample),
                   "SumOfCoherentIntens");
  }
  return tr;
}

void IncoherentIntensity::GetParameters(ComPWA::ParameterList &list) {
  list.AddParameter(_strength);
  for (auto i : GetIntensities()) {
    i->GetParameters(list);
  }
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
