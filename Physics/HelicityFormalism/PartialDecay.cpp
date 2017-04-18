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
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::vector<int> stringToVectInt(std::string str) {
  std::vector<int> result;
  std::istringstream iStr(str);
  std::vector<std::string> stringFrag{std::istream_iterator<std::string>{iStr},
                                      std::istream_iterator<std::string>{}};
  for (auto i : stringFrag) {
    result.push_back(std::stoi(i));
  }
  return result;
}

std::shared_ptr<PartialDecay>
PartialDecay::Factory(const boost::property_tree::ptree &pt) {

  LOG(trace) << "PartialDecay::Factory() |";
  auto obj = std::make_shared<PartialDecay>();
  obj->SetName(pt.get<std::string>("<xmlattr>.Name", "empty"));

  auto mag = ComPWA::DoubleParameterFactory(pt.get_child("Magnitude"));
  obj->SetMagnitudePar(std::make_shared<DoubleParameter>(mag));
  auto phase = ComPWA::DoubleParameterFactory(pt.get_child("Phase"));
  obj->SetPhasePar(std::make_shared<DoubleParameter>(phase));

  // Read subSystem definition
  std::vector<int> recoilState;
  auto recoil =
      pt.get_optional<std::string>("RecoilSystem.<xmlattr>.FinalState");
  if (recoil)
    recoilState = stringToVectInt(recoil.get());

  auto decayProducts = pt.get_child("DecayProducts");
  if (decayProducts.size() != 2)
    throw std::runtime_error(
        "PartialDecay::Factory() | Number of decay productes != 2.");

  auto itr = decayProducts.begin();
  auto finalStateA =
      stringToVectInt(itr->second.get<std::string>("<xmlattr>.FinalState"));
  auto finalStateB =
      stringToVectInt((++itr)->second.get<std::string>("<xmlattr>.FinalState"));
  SubSystem subSys(recoilState, finalStateA, finalStateB);

  // Create WignerD object
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
    dynObj = RelativisticBreitWigner::Factory(pt);
  } else {
    throw std::runtime_error("PartialDecay::Factory() | Unknown decay type " +
                             decayType + "!");
  }
  
  
  obj->SetDynamicalFunction(dynObj);

  //make sure dynamical function is created and set first
  obj->SetSubSystem(subSys);
  
  return obj;
}

/**! Setup function tree */
std::shared_ptr<FunctionTree> PartialDecay::GetTree(ParameterList &sample,
                                                    ParameterList &phspSample,
                                                    ParameterList &toySample,
                                                    std::string suffix) {
  std::shared_ptr<FunctionTree> tr(new FunctionTree());
  tr->createHead("Resonance(" + GetName() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)));
  tr->createNode("Strength",
                 std::shared_ptr<Strategy>(new Complexify(ParType::COMPLEX)),
                 "Resonance(" + GetName() + ")" + suffix);
  tr->createLeaf("Magnitude", _magnitude, "Strength");
  tr->createLeaf("Phase", _phase, "Strength");
  tr->createLeaf("PreFactor", _preFactor,
                 "Resonance(" + GetName() + ")" + suffix);

  tr->insertTree(_angD->GetTree(sample, (_dataPos*3)+1, (_dataPos*3)+2),
                 "Resonance(" + GetName() + ")" + suffix);
  tr->insertTree(_dynamic->GetTree(sample, toySample, (_dataPos*3)),
                 "Resonance(" + GetName() + ")" + suffix);

  return tr;
}
} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
