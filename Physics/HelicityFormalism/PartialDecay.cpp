//
//  PartialDecay.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 21/02/2017.
//
//

#include <sstream>

#include "Physics/DecayDynamics/RelativisticBreitWigner.hpp"
#include "Physics/DecayDynamics/AmpFlatteRes.hpp"
#include "Physics/DecayDynamics/NonResonant.hpp"

#include "Physics/HelicityFormalism/PartialDecay.hpp"
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
  obj->SetMagnitudeParameter(std::make_shared<DoubleParameter>(mag));
  auto phase = ComPWA::DoubleParameterFactory(pt.get_child("Phase"));
  obj->SetPhaseParameter(std::make_shared<DoubleParameter>(phase));

  // Read subSystem definition
  std::vector<int> recoilState;
  auto recoil =
      pt.get_optional<std::string>("RecoilSystem.<xmlattr>.FinalState");
  if (recoil){
    recoilState = stringToVectInt(recoil.get());
  }

  auto decayProducts = pt.get_child("DecayProducts");
  std::vector<std::vector<int>> finalStates;
  std::vector<std::string> finalStatesNames;
  for (auto i : decayProducts){
    finalStates.push_back(
        stringToVectInt(i.second.get<std::string>("<xmlattr>.FinalState")));
    finalStatesNames.push_back( i.second.get<std::string>("<xmlattr>.Name") );
  }
  SubSystem subSys(recoilState, finalStates);
  subSys.SetFinalStatesNames(finalStatesNames);

  auto dynObj = std::shared_ptr<DecayDynamics::AbstractDynamicalFunction>();
  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  
  if (finalStates.size() == 2) {

    // Create WignerD object
    obj->SetWignerD(
        ComPWA::Physics::HelicityFormalism::AmpWignerD::Factory(pt));

    auto partProp = PhysConst::Instance()->FindParticle(name);
    std::string decayType = partProp.GetDecayType();

    if (decayType == "stable") {
      throw std::runtime_error("PartialDecay::Factory() | Stable particle is "
                               "given as mother particle of a decay. Makes no "
                               "sense!");
    } else if (decayType == "relativisticBreitWigner") {
      dynObj = DecayDynamics::RelativisticBreitWigner::Factory(pt);
    } else if (decayType == "flatte") {
      dynObj = DecayDynamics::AmpFlatteRes::Factory(pt);
    } else {
      throw std::runtime_error("PartialDecay::Factory() | Unknown decay type " +
                               decayType + "!");
    }

    // make sure dynamical function is created and set first
  } else {
    dynObj = std::shared_ptr<DecayDynamics::AbstractDynamicalFunction>(new DecayDynamics::NonResonant);
    dynObj->SetName(name);
    // We assume the we have a multi-body decay and assume that the decay
    // proceeds via constant (non-resonant) dynamics
    obj->SetWignerD(std::shared_ptr<AmpWignerD>(new AmpWignerD()));
  }
  
  obj->SetDynamicalFunction(dynObj);
  obj->SetSubSystem(subSys);

  return std::static_pointer_cast<Resonance>(obj);
}

boost::property_tree::ptree PartialDecay::Save(std::shared_ptr<Resonance> res) {

  auto obj = std::static_pointer_cast<PartialDecay>(res);
  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->GetName());
  pt.add_child("Magnitude",
               ComPWA::DoubleParameterSave(*obj->GetMagnitudeParameter().get()));
  pt.add_child("Phase", ComPWA::DoubleParameterSave(*obj->GetPhaseParameter().get()));

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

  // Information daugher final state A
  auto finalS = obj->GetSubSystem().GetFinalStates();
  auto finalSNames = obj->GetSubSystem().GetFinalStatesNames();
  for (int j=0; j< finalS.size(); j++) {
    std::string strA;
    if (finalS.at(j).size()) {
      for (auto i = 0; i < finalS.at(j).size(); i++) {
        strA += std::to_string(finalS.at(j).at(i));
        if (i < finalS.at(j).size() - 1)
          strA += " ";
      }
    }

    chTrA.put("<xmlattr>.Name", finalSNames.at(j));
    chTrA.put("<xmlattr>.FinalState", strA);
    chTrA.put("<xmlattr>.Helicity",
              std::to_string((int)obj->GetWignerD()->GetHelicities().first));
    daughterTr.add_child("Particle", chTrA);
  }

  pt.add_child("DecayProducts", daughterTr);
  return pt;
}

/**! Setup function tree */
std::shared_ptr<FunctionTree> PartialDecay::GetTree(const ParameterList &sample,
                                                    const ParameterList &toySample,
                                                    std::string suffix) {

  int phspSampleSize = toySample.GetMultiDouble(0)->GetNValues();
  double phspVol = Kinematics::Instance()->GetPhspVolume();

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
  tr->insertTree(_angD->GetTree(sample, _dataPos + 1, _dataPos + 2),
                 "Resonance(" + GetName() + ")" + suffix);
  tr->insertTree(_dynamic->GetTree(sample,_dataPos),
                 "Resonance(" + GetName() + ")" + suffix);

  tr->recalculate();
  tr->createNode("Normalization",
                 std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)),
                 "Resonance(" + GetName() + ")" + suffix); // 1/normLH
  tr->createNode("SqrtIntegral",
                 std::shared_ptr<Strategy>(new SquareRoot(ParType::DOUBLE)),
                 "Normalization");
  tr->createNode("Integral",
                 std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
                 "SqrtIntegral");
  tr->createLeaf("PhspVolume", phspVol, "Integral");
  tr->createLeaf("InverseSampleSize", 1 / ((double)phspSampleSize), "Integral");
  tr->createNode("Sum", std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                 "Integral");
  tr->createNode("Intensity",
                 std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                 "Sum", phspSampleSize,
                 false); //|T_{ev}|^2
  tr->createNode("mult",
                 std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)),
                 "Intensity");
  tr->insertTree(_angD->GetTree(toySample, _dataPos + 1, _dataPos + 2, "_norm"),
                 "mult" + suffix);
  tr->insertTree(_dynamic->GetTree(toySample, _dataPos, "_norm"), "mult" + suffix);

  tr->recalculate();
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
