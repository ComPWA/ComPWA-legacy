// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <sstream>

#include "Physics/DecayDynamics/RelativisticBreitWigner.hpp"
#include "Physics/DecayDynamics/AmpFlatteRes.hpp"
#include "Physics/DecayDynamics/NonResonant.hpp"

#include "Physics/HelicityFormalism/PartialDecay.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::shared_ptr<Resonance>
PartialDecay::Factory(std::shared_ptr<PartList> partL,
                      std::shared_ptr<Kinematics> kin,
                      const boost::property_tree::ptree &pt) {

  LOG(trace) << "PartialDecay::Factory() |";
  SubSystem subSys = SubSystemFactory(pt);

  auto obj = std::make_shared<PartialDecay>(subSys);
  obj->SetName(pt.get<std::string>("<xmlattr>.Name", "empty"));

  std::shared_ptr<DoubleParameter> mag, phase;
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter") {
      if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude") {
        auto tmp = ComPWA::DoubleParameterFactory(v.second);
        mag = std::make_shared<DoubleParameter>(tmp);
      }
      if (v.second.get<std::string>("<xmlattr>.Type") == "Phase") {
        auto tmp = ComPWA::DoubleParameterFactory(v.second);
        phase = std::make_shared<DoubleParameter>(tmp);
      }
    } else {
      // ignored further settings. Should we throw an error?
    }
  }

  if (mag)
    obj->SetMagnitudeParameter(mag);
  else
    throw BadParameter(
        "PartialDecay::Factory() | No magnitude parameter found.");

  if (phase)
    obj->SetPhaseParameter(phase);
  else
    throw BadParameter(
        "PartialDecay::Factory() | No phase parameter found.");

  auto dynObj = std::shared_ptr<DecayDynamics::AbstractDynamicalFunction>();
  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");

  // Two-body decay
  if (subSys.GetFinalStates().size() == 2) {
    // Create WignerD object
    obj->SetWignerD(
        ComPWA::Physics::HelicityFormalism::AmpWignerD::Factory(partL, pt));

    auto partProp = partL->find(name)->second;
    std::string decayType = partProp.GetDecayType();

    if (decayType == "stable") {
      throw std::runtime_error("PartialDecay::Factory() | Stable particle is "
                               "given as mother particle of a decay. Makes no "
                               "sense!");
    } else if (decayType == "relativisticBreitWigner") {
      dynObj = DecayDynamics::RelativisticBreitWigner::Factory(partL, pt);
    } else if (decayType == "flatte") {
      dynObj = DecayDynamics::AmpFlatteRes::Factory(partL, pt);
    } else {
      throw std::runtime_error("PartialDecay::Factory() | Unknown decay type " +
                               decayType + "!");
    }

    // make sure dynamical function is created and set first
  } else { // Multi-body decay
    dynObj = std::shared_ptr<DecayDynamics::AbstractDynamicalFunction>(
        new DecayDynamics::NonResonant);
    dynObj->SetName(name);
    // We assume the we have a multi-body decay and assume that the decay
    // proceeds via constant (non-resonant) dynamics
    obj->SetWignerD(std::shared_ptr<AmpWignerD>(new AmpWignerD()));
  }

  obj->SetDynamicalFunction(dynObj);
  obj->SetDataPosition(
      3 *
      std::dynamic_pointer_cast<HelicityKinematics>(kin)->GetDataID(subSys));

  obj->SetPhspVolume(kin->GetPhspVolume());

  return std::static_pointer_cast<Resonance>(obj);
}

boost::property_tree::ptree PartialDecay::Save(std::shared_ptr<Resonance> res) {

  auto obj = std::static_pointer_cast<PartialDecay>(res);
  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->GetName());
  pt.add_child("Magnitude", ComPWA::DoubleParameterSave(
                                *obj->GetMagnitudeParameter().get()));
  pt.add_child("Phase",
               ComPWA::DoubleParameterSave(*obj->GetPhaseParameter().get()));

  pt.put("DecayParticle.<xmlattr>.Name",
         obj->GetDynamicalFunction()->GetName());
  pt.put("DecayParticle.<xmlattr>.Helicity", obj->GetWignerD()->GetMu());

  auto subTr = SubSystemSave(obj->GetSubSystem());
  pt.add_child("SubSystem", subTr);

  return pt;
}

bool PartialDecay::CheckModified() const {
  if (Resonance::CheckModified())
    return true;
  if (_dynamic->CheckModified()) {
    const_cast<double &>(_current_integral) = Integral();
    _dynamic->SetModified(false);
    return true;
  }
  return false;
}

double PartialDecay::GetNormalization() const {
  if (_dynamic->CheckModified() || !_current_integral)
    const_cast<double &>(_current_integral) = Integral();
  _dynamic->SetModified(false);
  assert(_current_integral != 0.0);
  return 1 / std::sqrt(_current_integral);
}

std::shared_ptr<FunctionTree>
PartialDecay::GetTree(std::shared_ptr<Kinematics> kin,
                      const ParameterList &sample,
                      const ParameterList &toySample, std::string suffix) {

  int phspSampleSize = toySample.GetMultiDouble(0)->GetNValues();

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
  tr->insertTree(_dynamic->GetTree(sample, _dataPos),
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
  tr->createLeaf("PhspVolume", phspVolume_, "Integral");
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
  tr->insertTree(_dynamic->GetTree(toySample, _dataPos, "_norm"),
                 "mult" + suffix);

  tr->recalculate();
  return tr;
}

void PartialDecay::GetParameters(ParameterList &list) {
  Resonance::GetParameters(list);
  //    _angD->GetParameters(list);
  _dynamic->GetParameters(list);
}

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA
