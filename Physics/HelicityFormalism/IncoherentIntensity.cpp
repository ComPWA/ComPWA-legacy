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
  obj->SetName(
      pt.get<std::string>("IncoherentIntensity.<xmlattr>.Name", "empty"));

  for (const auto &v : pt.get_child("IncoherentIntensity")) {
    if (v.first == "CoherentIntensity")
      obj->Add(ComPWA::Physics::HelicityFormalism::CoherentIntensity::Factory(
          v.second));
  }
  return obj;
}

std::shared_ptr<ComPWA::FunctionTree> IncoherentIntensity::GetTree(
    ComPWA::ParameterList &sample, ComPWA::ParameterList &phspSample,
    ComPWA::ParameterList &toySample, std::string suffix) {

  std::shared_ptr<FunctionTree> tr(new FunctionTree());
  
  tr->createHead("IncoherentSum",
                 std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)));
  tr->createLeaf("Strength", _strength, "IncoherentSum");
  tr->createNode("SumCoherent",
                 std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),"IncoherentSum");
  for (auto i : _intens) {
    tr->insertTree(i->GetTree(sample, phspSample, toySample), "SumCoherent");
  }
  return tr;
}
  
} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
