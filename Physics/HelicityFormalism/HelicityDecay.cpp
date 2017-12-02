// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <sstream>

#include "Physics/DecayDynamics/RelativisticBreitWigner.hpp"
#include "Physics/DecayDynamics/AmpFlatteRes.hpp"
#include "Physics/DecayDynamics/NonResonant.hpp"

#include "Physics/HelicityFormalism/HelicityDecay.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::shared_ptr<PartialAmplitude>
HelicityDecay::Factory(std::shared_ptr<PartList> partL,
                      std::shared_ptr<Kinematics> kin,
                      const boost::property_tree::ptree &pt) {

  LOG(trace) << "HelicityDecay::Factory() |";
  SubSystem subSys = SubSystemFactory(pt);

  // Create object. Setting dataPos to invalid, will be set later
  auto obj = std::make_shared<HelicityDecay>(-1, subSys);
  obj->setName(pt.get<std::string>("<xmlattr>.Name", "empty"));

  std::shared_ptr<DoubleParameter> mag, phase;
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter") {
      if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude") {
        auto tmp = DoubleParameter();
        tmp.load(v.second);
        mag = std::make_shared<DoubleParameter>(tmp);
      }
      if (v.second.get<std::string>("<xmlattr>.Type") == "Phase") {
        auto tmp = DoubleParameter();
        tmp.load(v.second);
        phase = std::make_shared<DoubleParameter>(tmp);
      }
    } else {
      // ignored further settings. Should we throw an error?
    }
  }

  if (mag)
    obj->setMagnitudeParameter(mag);
  else {
    obj->setMagnitudeParameter(std::make_shared<ComPWA::DoubleParameter>("", 1.0));
    obj->magnitudeParameter()->fixParameter(true);
  }
  
  if (phase)
    obj->setPhaseParameter(phase);
  else {
    obj->setPhaseParameter(std::make_shared<ComPWA::DoubleParameter>("", 0.0));
    obj->phaseParameter()->fixParameter(true);
  }

  auto dynObj = std::shared_ptr<DecayDynamics::AbstractDynamicalFunction>();
  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");

  // Two-body decay
  if (subSys.GetFinalStates().size() == 2) {
    // Create WignerD object
    obj->setWignerD(
        ComPWA::Physics::HelicityFormalism::AmpWignerD::Factory(partL, pt));

    auto partProp = partL->find(name)->second;
    std::string decayType = partProp.GetDecayType();

    if (decayType == "stable") {
      throw std::runtime_error("HelicityDecay::Factory() | Stable particle is "
                               "given as mother particle of a decay. Makes no "
                               "sense!");
    } else if (decayType == "relativisticBreitWigner") {
      dynObj = DecayDynamics::RelativisticBreitWigner::Factory(partL, pt);
    } else if (decayType == "flatte") {
      dynObj = DecayDynamics::AmpFlatteRes::Factory(partL, pt);
    } else {
      throw std::runtime_error("HelicityDecay::Factory() | Unknown decay type " +
                               decayType + "!");
    }

    // make sure dynamical function is created and set first
  } else { // Multi-body decay
    dynObj = std::shared_ptr<DecayDynamics::AbstractDynamicalFunction>(
        new DecayDynamics::NonResonant);
    dynObj->setName(name);
    // We assume the we have a multi-body decay and assume that the decay
    // proceeds via constant (non-resonant) dynamics
    obj->setWignerD(std::shared_ptr<AmpWignerD>(new AmpWignerD()));
  }

  obj->setDynamicalFunction(dynObj);
  obj->setDataPosition(
      3 *
      std::dynamic_pointer_cast<HelicityKinematics>(kin)->dataID(subSys));

  obj->setPhspVolume(kin->phspVolume());

  return std::static_pointer_cast<PartialAmplitude>(obj);
}

boost::property_tree::ptree HelicityDecay::Save(std::shared_ptr<PartialAmplitude> res) {

  auto obj = std::static_pointer_cast<HelicityDecay>(res);
  auto pt = SubSystemSave(obj->subSystem());
  pt.put<std::string>("<xmlattr>.Name", obj->name());
  
  boost::property_tree::ptree tmp = obj->magnitudeParameter()->save();
  tmp.put("<xmlattr>.Type", "Magnitude");
  tmp.put("<xmlattr>.Class", "Double");
  pt.add_child("Parameter", tmp);
  
  tmp = obj->phaseParameter()->save();
  tmp.put("<xmlattr>.Type", "Phase");
  tmp.put("<xmlattr>.Class", "Double");
  pt.add_child("Parameter", tmp);
  
  pt.put("DecayParticle.<xmlattr>.Name",
         obj->dynamicalFunction()->name());
  pt.put("DecayParticle.<xmlattr>.Helicity", obj->wignerD()->GetMu());

  //TODO: put helicities of daughter particles
  return pt;
}

bool HelicityDecay::isModified() const {
  if (PartialAmplitude::isModified())
    return true;
  if (DynamicFcn->CheckModified()) {
    const_cast<double &>(CurrentIntegral) = integral();
    DynamicFcn->SetModified(false);
    return true;
  }
  return false;
}

double HelicityDecay::normalization() const {
  if (DynamicFcn->CheckModified() || !CurrentIntegral)
    const_cast<double &>(CurrentIntegral) = integral();
  DynamicFcn->SetModified(false);
  assert(CurrentIntegral != 0.0);
  return 1 / std::sqrt(CurrentIntegral);
}

std::shared_ptr<FunctionTree>
HelicityDecay::tree(std::shared_ptr<Kinematics> kin,
                      const ParameterList &sample,
                      const ParameterList &toySample, std::string suffix) {

  int phspSampleSize = toySample.GetMultiDouble(0)->numValues();

  std::shared_ptr<FunctionTree> tr(new FunctionTree());
  tr->createHead("PartialAmplitude(" + name() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)));
  tr->createNode("Strength",
                 std::shared_ptr<Strategy>(new Complexify(ParType::COMPLEX)),
                 "PartialAmplitude(" + name() + ")" + suffix);
  tr->createLeaf("Magnitude", _magnitude, "Strength");
  tr->createLeaf("Phase", _phase, "Strength");
  tr->createLeaf("PreFactor", _preFactor,
                 "PartialAmplitude(" + name() + ")" + suffix);
  tr->insertTree(AngularDist->GetTree(sample, DataPosition + 1, DataPosition + 2),
                 "PartialAmplitude(" + name() + ")" + suffix);
  tr->insertTree(DynamicFcn->GetTree(sample, DataPosition),
                 "PartialAmplitude(" + name() + ")" + suffix);

  tr->recalculate();
  tr->createNode("Normalization",
                 std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)),
                 "PartialAmplitude(" + name() + ")" + suffix); // 1/normLH
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
  tr->insertTree(AngularDist->GetTree(toySample, DataPosition + 1, DataPosition + 2, "_norm"),
                 "mult" + suffix);
  tr->insertTree(DynamicFcn->GetTree(toySample, DataPosition, "_norm"),
                 "mult" + suffix);

  tr->recalculate();
  return tr;
}

void HelicityDecay::parameters(ParameterList &list) {
  PartialAmplitude::parameters(list);
  //    AngularDist->GetParameters(list);
  DynamicFcn->GetParameters(list);
}

void HelicityDecay::updateParameters(const ParameterList &list){

  // Try to update magnitude
  std::shared_ptr<DoubleParameter> mag;
  try {
    mag = list.GetDoubleParameter(_magnitude->name());
  } catch (std::exception &ex) {
  }
  if (mag)
    _magnitude->updateParameter(mag);
  std::shared_ptr<DoubleParameter> phase;
  
  // Try to update phase
  try {
    phase = list.GetDoubleParameter(_phase->name());
  } catch (std::exception &ex) {
  }
  if (phase)
    _phase->updateParameter(phase);

  DynamicFcn->updateParameters(list);
  
  return;

}

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA
