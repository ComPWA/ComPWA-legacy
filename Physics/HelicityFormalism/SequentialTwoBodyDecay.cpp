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
    obj->SetName( pt.get<std::string>("Amplitude.<xmlattr>.Name","empty") );
    
      auto mag= ComPWA::DoubleParameterFactory( pt.get_child("Magnitude") );
      obj->SetMagnitude( std::make_shared<DoubleParameter>(mag) );
      auto phase= ComPWA::DoubleParameterFactory( pt.get_child("Phase") );
      obj->SetPhase( std::make_shared<DoubleParameter>(phase) );
    
    for(const auto& v : pt.get_child("") ){
      if( v.first == "Resonance" )
        obj->Add( PartialDecay::Factory(v.second) );
    }
    return obj;
  }
} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
