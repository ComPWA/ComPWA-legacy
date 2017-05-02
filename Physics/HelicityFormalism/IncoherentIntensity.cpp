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
      pt.get<std::string>("<xmlattr>.Name", "empty"));

  auto ptCh = pt.get_child_optional("Strength");
  if (ptCh) {
    auto strength = ComPWA::DoubleParameterFactory(ptCh.get());
    obj->SetStrengthParameter(std::make_shared<DoubleParameter>(strength));
  } else {
    obj->SetStrengthParameter(std::make_shared<ComPWA::DoubleParameter>("", 1.0));
    obj->GetStrengthParameter()->SetParameterFixed();
  }
  
  for (const auto &v : pt.get_child("")) {
    if (v.first == "CoherentIntensity")
      obj->Add(ComPWA::Physics::HelicityFormalism::CoherentIntensity::Factory(
          v.second));
  }
  return obj;
}

boost::property_tree::ptree
IncoherentIntensity::Save(std::shared_ptr<IncoherentIntensity> obj) {

  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name",obj->GetName());
  pt.add_child("Strength", ComPWA::DoubleParameterSave(*obj->GetStrengthParameter().get()));
  for( auto i : obj->GetIntensities() ) {
    pt.add_child("CoherentIntensity", CoherentIntensity::Save(i));
  }
  return pt;
}

std::shared_ptr<ComPWA::FunctionTree> IncoherentIntensity::GetTree(
    const ComPWA::ParameterList &sample, const ComPWA::ParameterList &phspSample,
    const ComPWA::ParameterList &toySample, std::string suffix) {

  std::shared_ptr<FunctionTree> tr(new FunctionTree());

  tr->createHead("IncoherentIntens(" + GetName() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  tr->createLeaf("Strength", _strength,
                 "IncoherentIntens(" + GetName() + ")" + suffix);
  tr->createNode("SumOfCoherentIntens",
                 std::shared_ptr<Strategy>(new AddAll(ParType::MDOUBLE)),
                 "IncoherentIntens(" + GetName() + ")" + suffix);
  for (auto i : _intens) {
    tr->insertTree(i->GetTree(sample, phspSample, toySample),
                   "SumOfCoherentIntens");
  }
  return tr;
}
  
  void IncoherentIntensity::GetParameters(ComPWA::ParameterList &list){
    
    list.AddParameter(GetStrengthParameter());
    for( auto i : GetIntensities() ) {
      i->GetParameters(list);
    }
    
  }
  
} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
