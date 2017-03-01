//
//  IncoherentIntensity.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 21/02/2017.
//
//

#include <stdio.h>
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::shared_ptr<IncoherentIntensity>
IncoherentIntensity::Factory(const boost::property_tree::ptree &pt) {

  LOG(trace) << " IncoherentIntensity::Factory() | Construction....";
  auto obj = std::make_shared<IncoherentIntensity>();
  obj->SetName(pt.get<std::string>("AmpIntensity.<xmlattr>.Name", "empty"));

  for (const auto &v : pt.get_child("CoherentIntensity")) {
    obj->Add(ComPWA::Physics::HelicityFormalism::CoherentIntensity::Factory(
        v.second));
  }
  return obj;
}
  
}
}
}
