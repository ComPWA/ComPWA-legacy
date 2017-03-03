//
//  PartialDecay.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 21/02/2017.
//
//

#include <stdio.h>

#include "Physics/HelicityFormalism/PartialDecay.hpp"
#include "Physics/HelicityFormalism/RelativisticBreitWigner.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::shared_ptr<PartialDecay>
PartialDecay::Factory(const boost::property_tree::ptree &pt) {

  LOG(trace) << "PartialDecay::Factory() |";
  auto obj = std::make_shared<PartialDecay>();
  obj->SetName(pt.get<string>("<xmlattr>.Name", "empty"));

  auto mag = ComPWA::DoubleParameterFactory(pt.get_child("Magnitude"));
  obj->SetMagnitudePar(std::make_shared<DoubleParameter>(mag));
  auto phase = ComPWA::DoubleParameterFactory(pt.get_child("Phase"));
  obj->SetPhasePar(std::make_shared<DoubleParameter>(phase));

  obj->SetWignerD(ComPWA::Physics::HelicityFormalism::AmpWignerD::Factory(pt));

  auto dynObj = std::shared_ptr<AbstractDynamicalFunction>();
  int id = pt.get<double>("DecayParticle.<xmlattr>.Id");
  auto partProp = PhysConst::Instance()->FindParticle(id);
  std::string decayType = partProp.GetDecayType();
  
  
  
  if (decayType == "stable") {
    throw std::runtime_error("PartialDecay::Factory() | Stable particle is "
                             "given as mother particle of a decay. Makes no "
                             "sense!");
  } else if (decayType == "relativisticBreitWigner") {
    dynObj = RelativisticBreitWigner::Factory( pt );
  } else {
    throw std::runtime_error("PartialDecay::Factory() | Unknown decay type " +
                             decayType + "!");
  }
  obj->SetDynamicalFunction(dynObj);

  return obj;
}
}
}
}
