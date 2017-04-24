  //
//  SequentialTwoBodyDecay.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 01/03/2017.
//
//

#include <stdio.h>
#include "Physics/HelicityFormalism/SequentialTwoBodyDecay.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::shared_ptr<SequentialTwoBodyDecay>
SequentialTwoBodyDecay::Factory(const boost::property_tree::ptree &pt) {
  LOG(trace) << " SequentialTwoBodyDecay::Factory() | Construction....";
  auto obj = std::make_shared<SequentialTwoBodyDecay>();
  obj->SetName(pt.get<std::string>("<xmlattr>.Name", "empty"));

  auto mag = ComPWA::DoubleParameterFactory(pt.get_child("Magnitude"));
  obj->SetMagnitude(std::make_shared<DoubleParameter>(mag));
  auto phase = ComPWA::DoubleParameterFactory(pt.get_child("Phase"));
  obj->SetPhase(std::make_shared<DoubleParameter>(phase));

  for (const auto &v : pt.get_child("")) {
    if (v.first == "Resonance")
      obj->Add(PartialDecay::Factory(v.second));
  }
  return obj;
}
  
boost::property_tree::ptree
SequentialTwoBodyDecay::Save(std::shared_ptr<SequentialTwoBodyDecay> obj) {

  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name",obj->GetName());
  pt.add_child("Magnitude",
               ComPWA::DoubleParameterSave(*obj->GetMagnitude().get()));
  pt.add_child("Phase", ComPWA::DoubleParameterSave(*obj->GetPhase().get()));
  
  for( auto i : obj->GetDecays() ) {
    pt.add_child("Resonance", PartialDecay::Save(i));
  }
  return pt;
}
  
/**! Setup function tree */
std::shared_ptr<FunctionTree>
SequentialTwoBodyDecay::GetTree(ParameterList &sample,
                                ParameterList &phspSample,
                                ParameterList &toySample, std::string suffix) {

  std::shared_ptr<FunctionTree> tr(new FunctionTree());
  tr->createHead("Amplitude("+GetName()+")"+suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)));

  for (auto i : _partDecays) {
    std::shared_ptr<FunctionTree> resTree =
        i->GetTree(sample, phspSample, phspSample, "");
    if (!resTree->sanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Amplitude tree didn't pass sanity check!");
    resTree->recalculate();
    tr->insertTree(resTree, "Amplitude("+GetName()+")"+suffix);
  }

  return tr;
};
  
void SequentialTwoBodyDecay::GetParameters(ParameterList &list) {
  Amplitude::GetParameters(list);
  for( auto i: _partDecays){
    i->GetParameters(list);
  }
}
  
} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
