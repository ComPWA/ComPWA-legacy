// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.
//

#include <numeric>

#include "EvtGenIF.hpp"
#include "Physics/BuilderXML.hpp"

namespace ComPWA {
namespace Physics {
namespace EvtGen {

using ComPWA::Physics::SubSystem;

void EvtGenIF::addResonance(const std::string &name, double m0, double g0,
                            double spin, const SubSystem &subsys) {
  // LOG(debug) << "EvtGenIF::addResonance starts. Name: " << name << " mass: "
  // << m0 << " width: " << g0 << " spin: " << spin;
  EvtCyclic3::Pair pairAng;
  EvtCyclic3::Pair pairRes;
  LOG(DEBUG) << "EvtGenIF::addResonance num finalstate00: "
             << subsys.getFinalStates().at(0).at(0);
  LOG(DEBUG) << "EvtGenIF::addResonance num finalstate01: "
             << subsys.getFinalStates().at(1).at(0);
  switch (subsys.getFinalStates().at(0).at(0) +
          subsys.getFinalStates().at(1).at(0)) {
  case 3: // 12:
    pairAng = EvtCyclic3::Pair::AB;
    pairRes = EvtCyclic3::Pair::BC;
    break;
  case 4: // 13:
    pairAng = EvtCyclic3::Pair::AC;
    pairRes = EvtCyclic3::Pair::AB;
    break;
  default: // 23
    pairAng = EvtCyclic3::Pair::BC;
    pairRes = EvtCyclic3::Pair::AC;
  }
  auto evtspin = EvtSpinType::spintype::SCALAR;
  if (spin == 1)
    evtspin = EvtSpinType::spintype::VECTOR;
  else if (spin == 2)
    evtspin = EvtSpinType::spintype::TENSOR;
  else
    throw std::runtime_error(
        "EvtGenIF::addResonance() | Spins higher than two are not supported.");

  LOG(DEBUG) << "EvtGenIF::addResonance: " << name;
  LOG(DEBUG) << "EvtGenIF::addResonance mass: " << m0;
  LOG(DEBUG) << "EvtGenIF::addResonance width: " << g0;
  LOG(DEBUG) << "EvtGenIF::addResonance spin: " << spin;
  EvtDalitzReso reso(name, DalitzPlot, pairAng, pairRes, evtspin, m0, g0,
                     EvtDalitzReso::NumType::RBW_ZEMACH);
  Resos.push_back(reso);

  evtPars[std::string(name + "_mass")] =
      std::shared_ptr<ComPWA::FunctionTree::FitParameter>(
          new ComPWA::FunctionTree::FitParameter(name + "_mass", m0));
  evtPars.at(std::string(name + "_mass"))->fixParameter(false);
  evtPars[std::string(name + "_width")] =
      std::shared_ptr<ComPWA::FunctionTree::FitParameter>(
          new ComPWA::FunctionTree::FitParameter(name + "_width", g0));
  evtPars.at(std::string(name + "_width"))->fixParameter(true);
  // LOG(debug) << "EvtGenIF::addResonance finishes";
}

void EvtGenIF::addHeliResonance(const boost::property_tree::ptree &pt,
                                const ComPWA::ParticleList &partL) {
  // LOG(debug) << "EvtGenIF::addHeliResonance starts";
  std::string resoName = pt.get<std::string>("<xmlattr>.Name", "empty");
  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  auto parti = findParticle(partL, name);

  double J = parti.getQuantumNumber<double>("Spin");

  auto DecayInfo = ComPWA::Physics::extractDecayInfo(pt);

  // LOG(debug) << "EvtGenIF::addHeliResonance two-body decay";
  // Two-body decay
  if (DecayInfo.SubSys.getFinalStates().size() == 2) {
    // Create WignerD object
    // AngularDist = std::make_shared<HelicityFormalism::AmpWignerD>(
    //     J, mu, DecayHelicities.first - DecayHelicities.second);

    auto partProp = findParticle(partL, name);
    std::string decayType = partProp.getDecayType();

    if (decayType == "stable") {
      throw std::runtime_error(
          "EvtGenIF::addHeliResonance() | Stable particle is "
          "given as mother particle of a decay. Makes no "
          "sense!");
    } else if (decayType == "relativisticBreitWigner") {
      LOG(DEBUG) << "EvtGenIF::addHeliResonance add resonance";
      double width = 0.1;
      for (const auto &v : partProp.getDecayInfo().get_child("")) {
        if (v.first != "Parameter")
          continue;
        std::string type = v.second.get<std::string>("<xmlattr>.Type");
        if (type == "Width") {
          // SetWidthParameter(std::make_shared<FitParameter>(v.second));
          width = v.second.get<double>("Value"); // TODO: funzt nicht
        } else if (type == "MesonRadius") {
          // TODO
        }
      }
      addResonance(name, parti.getMass().Value, width, J, DecayInfo.SubSys);
      // DynamicFcn = std::make_shared<DecayDynamics::RelativisticBreitWigner>(
      //  name, DecayProducts, partL);
    } else if (decayType == "flatte") {
      // TODO
      // DynamicFcn = std::make_shared<DecayDynamics::AmpFlatteRes>(
      //  name, DecayProducts, partL);
    } else {
      throw std::runtime_error(
          "EvtGenIF::addHeliResonance() | Unknown decay type " + decayType +
          "!");
    }
    LOG(DEBUG) << "EvtGenIF::addHeliResonance finished";
    // make sure dynamical function is created and set first
  }
}

void EvtGenIF::addResonances(const boost::property_tree::ptree &pt,
                             std::shared_ptr<DalitzKinematics> kin,
                             const ComPWA::ParticleList &partL) {
  if (pt.get<std::string>("<xmlattr>.Class") != "Incoherent")
    throw BadConfig("IncoherentIntensity::Factory() | Property tree seems to "
                    "not containt a configuration for an "
                    "IncoherentIntensity!");

  setPhspVolume(kin->phspVolume());

  LOG(DEBUG) << "EvtGenIF::addResoances starts";

  for (const auto &v : pt.get_child("")) {
    if (v.first == "Intensity" &&
        v.second.get<std::string>("<xmlattr>.Class") == "Coherent") {

      for (const auto &w : v.second.get_child("")) {
        if (w.first == "Amplitude" &&
            w.second.get<std::string>("<xmlattr>.Class") ==
                "SequentialPartialAmplitude") {
          // addAmplitude(std::make_shared<SequentialPartialAmplitude>(partL,
          // kin, v.second));
          std::string name = pt.get<std::string>("<xmlattr>.Name", "empty");

          std::complex<double> preFactor = std::complex<double>(1, 0);
          for (const auto &x : w.second.get_child("")) {
            if (x.first == "Parameter") {

            } else if (x.first == "PartialAmplitude" &&
                       x.second.get<std::string>("<xmlattr>.Class") ==
                           "HelicityDecay") {
              LOG(DEBUG) << "EvtGenIF::addResoances add Heli-Resonance";
              addHeliResonance(x.second, partL);

            } else if (x.first == "PreFactor") {
              LOG(DEBUG) << "EvtGenIF::addResoances add factor";
              double r = x.second.get<double>("<xmlattr>.Magnitude");
              double p = x.second.get<double>("<xmlattr>.Phase");
              preFactor = std::polar(r, p);
            }
          }
        }
      }
    }
  }

  LOG(DEBUG) << "EvtGenIF::addResoances finished";
}

std::vector<double> EvtGenIF::evaluate(const ComPWA::DataMap &data) noexcept {
  std::vector<double> Results;
  for (size_t EventIndex = 0; EventIndex < data.at("mA").size(); ++EventIndex) {
    EvtDalitzPoint pnt(data.at("mA")[EventIndex], data.at("mB")[EventIndex],
                       data.at("mC")[EventIndex], data.at("qAB")[EventIndex],
                       data.at("qBC")[EventIndex], data.at("qCA")[EventIndex]);
    double result = 0;

    for (unsigned int i = 0; i < Resos.size(); ++i) {
      EvtDalitzReso tmp = Resos.at(i);
      std::string name = tmp.get_Name();

      tmp.set_Mass(evtPars.at(std::string(name + "_mass"))->value());
      tmp.set_Gamma(evtPars.at(std::string(name + "_width"))->value());

      result += abs2(tmp.evaluate(pnt)); // * normValues.at(i);
    }

    assert(!std::isnan(result) &&
           "IncoherentIntensity::Intensity() | Result is NaN!");
    assert(!std::isinf(result) &&
           "IncoherentIntensity::Intensity() | Result is inf!");

    Results.push_back(result);
  }
  return Results;
}

void EvtGenIF::updateParametersFrom(const std::vector<double> &Parameters) {
  size_t counter(0);
  for (auto i : evtPars) {
    if (!i.second->isFixed())
      i.second->setValue(Parameters[counter]);
    ++counter;
  }
}

std::vector<ComPWA::Parameter> EvtGenIF::getParameters() const {
  std::vector<ComPWA::Parameter> pars;
  for (auto i : evtPars) {
    pars.push_back(ComPWA::Parameter{i.second->name(), i.second->value()});
  }
  return pars;
}

} // namespace EvtGen
} // namespace Physics
} // namespace ComPWA
