//
//  PartialDecay.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 21/02/2017.
//
//

//#include <stdio.h>
#include <sstream>


#include "Physics/HelicityFormalism/PartialDecay.hpp"
#include "Physics/HelicityFormalism/RelativisticBreitWigner.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {
  
  std::vector<int> stringToVectInt( std::string str ){
    std::vector<int> result;
    std::istringstream iStr(str);
    vector<string> stringFrag{istream_iterator<string>{iStr},
      istream_iterator<string>{}};
    for( auto i: stringFrag ){
      result.push_back(std::stoi(i));
    }
    return result;
  }

std::shared_ptr<PartialDecay>
PartialDecay::Factory(const boost::property_tree::ptree &pt) {

  LOG(trace) << "PartialDecay::Factory() |";
  auto obj = std::make_shared<PartialDecay>();
  obj->SetName(pt.get<string>("<xmlattr>.Name", "empty"));

  
  auto mag = ComPWA::DoubleParameterFactory(pt.get_child("Magnitude"));
  obj->SetMagnitudePar(std::make_shared<DoubleParameter>(mag));
  auto phase = ComPWA::DoubleParameterFactory(pt.get_child("Phase"));
  obj->SetPhasePar(std::make_shared<DoubleParameter>(phase));

  std::vector<int> recoilState;
  auto recoil = pt.get_optional<std::string>("RecoilSystem.<xmlattr>.FinalState");
  if( recoil )
    recoilState = stringToVectInt(recoil.get());
  
  std::vector<int> finalState;
  auto decayProducts = pt.get_child("DecayProducts");
  for (auto i : decayProducts) {
    auto strFS = i.second.get<std::string>("<xmlattr>.FinalState");
    auto intFS = stringToVectInt(strFS);
    finalState.insert(finalState.end(), intFS.begin(), intFS.end());
  }
//  std::cout<<obj->GetName() <<" recoilState: ";
//  for( auto i : recoilState )
//    std::cout<<i <<" ";
//  std::cout<<std::endl<<" finalState: ";
//  for( auto i : finalState )
//    std::cout<<i <<" ";
//  std::cout<<std::endl;
  
  obj->SetWignerD(ComPWA::Physics::HelicityFormalism::AmpWignerD::Factory(pt));

  auto dynObj = std::shared_ptr<AbstractDynamicalFunction>();
  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  auto partProp = PhysConst::Instance()->FindParticle(name);
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

  //TODO: Need to get data poistion from Kinematics
//  int pos = Kinematics::Instance()
//  obj->SetDataPosition(mPos, cosThetaPos, phiPos)
  
  return obj;
}
  
}
}
}
