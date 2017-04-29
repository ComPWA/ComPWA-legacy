//
//  PartialDecay.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 21/02/2017.
//
//

#include <sstream>

#include "Physics/HelicityFormalism/PartialDecay.hpp"
#include "Physics/HelicityFormalism/RelativisticBreitWigner.hpp"
#include "Physics/HelicityFormalism/AmpFlatteRes.hpp"
#include "Physics/HelicityFormalism/NonResonant.hpp"
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

std::shared_ptr<Resonance>
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
  std::vector<std::vector<int>> finalStates;
  for (auto i : decayProducts)
    finalStates.push_back(
        stringToVectInt(i.second.get<std::string>("<xmlattr>.FinalState")));

  if (finalStates.size() == 2) {
    SubSystem subSys(recoilState, finalStates.at(0), finalStates.at(1));

    // Create WignerD object
    obj->SetWignerD(
        ComPWA::Physics::HelicityFormalism::AmpWignerD::Factory(pt));

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
    } else if (decayType == "flatte") {
      dynObj = AmpFlatteRes::Factory(pt);
    } else {
      throw std::runtime_error("PartialDecay::Factory() | Unknown decay type " +
                               decayType + "!");
    }
    dynObj->SetSubSystem(subSys);
    obj->SetDynamicalFunction(dynObj);

    // make sure dynamical function is created and set first
    obj->SetSubSystem(subSys);
  } else {
    // We assume the we have a multi-body decay and assume that the decay
    // proceeds via constant (non-resonant) dynamics
    obj->SetDynamicalFunction(std::shared_ptr<AbstractDynamicalFunction>(new NonResonant));
  }

  return std::static_pointer_cast<Resonance>(obj);
}

boost::property_tree::ptree PartialDecay::Save(std::shared_ptr<Resonance> res) {

  auto obj = std::static_pointer_cast<PartialDecay>(res);
  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->GetName());
  pt.add_child("Magnitude",
               ComPWA::DoubleParameterSave(*obj->GetMagnitudePar().get()));
  pt.add_child("Phase", ComPWA::DoubleParameterSave(*obj->GetPhasePar().get()));

  pt.put("DecayParticle.<xmlattr>.Name",
         obj->GetDynamicalFunction()->GetName());
  pt.put("DecayParticle.<xmlattr>.Helicity", obj->GetWignerD()->GetMu());

  auto recoilV = obj->GetSubSystem().GetRecoilState();
  if (recoilV.size()) {
    std::string recoilStr;
    for (auto i = 0; i < recoilV.size(); i++) {
      recoilStr += std::to_string(recoilV.at(i));
      if (i < recoilV.size() - 1)
        recoilStr += " ";
    }
    pt.put("RecoilSystem.<xmlattr>.FinalState", recoilStr);
  }

  boost::property_tree::ptree daughterTr, chTrA, chTrB;

  auto names = obj->GetDynamicalFunction()->GetDecayNames();

  // Information daugher final state A
  auto finalA = obj->GetSubSystem().GetFinalStateA();
  std::string strA;
  if (finalA.size()) {
    for (auto i = 0; i < finalA.size(); i++) {
      strA += std::to_string(finalA.at(i));
      if (i < finalA.size() - 1)
        strA += " ";
    }
  }

  chTrA.put("<xmlattr>.Name", names.first);
  chTrA.put("<xmlattr>.FinalState", strA);
  chTrA.put("<xmlattr>.Helicity",
            std::to_string((int)obj->GetWignerD()->GetHelicities().first));
  daughterTr.add_child("Particle", chTrA);

  // Information daugher final state B
  auto finalB = obj->GetSubSystem().GetFinalStateB();
  std::string strB;
  if (finalB.size()) {
    for (auto i = 0; i < finalB.size(); i++) {
      strB += std::to_string(finalB.at(i));
      if (i < finalB.size() - 1)
        strB += " ";
    }
  }

  chTrB.put("<xmlattr>.Name", names.second);
  chTrB.put("<xmlattr>.FinalState", strB);
  chTrB.put("<xmlattr>.Helicity",
            std::to_string((int)obj->GetWignerD()->GetHelicities().second));
  daughterTr.add_child("Particle", chTrB);

  pt.add_child("DecayProducts", daughterTr);
  return pt;
}

/**! Setup function tree */
std::shared_ptr<FunctionTree> PartialDecay::GetTree(ParameterList &sample,
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

  tr->insertTree(_angD->GetTree(sample, (_dataPos * 3) + 1, (_dataPos * 3) + 2),
                 "Resonance(" + GetName() + ")" + suffix);
  tr->insertTree(_dynamic->GetTree(sample, toySample),
                 "Resonance(" + GetName() + ")" + suffix);

  return tr;
}

void PartialDecay::GetParameters(ParameterList &list) {
  Resonance::GetParameters(list);
  //    _angD->GetParameters(list);
  _dynamic->GetParameters(list);
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
