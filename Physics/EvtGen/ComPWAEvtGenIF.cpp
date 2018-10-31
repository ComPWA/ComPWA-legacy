// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.
//

#include "Physics/EvtGen/ComPWAEvtGenIF.hpp"
#include <numeric>

using namespace ComPWA::Physics::EvtGenIF;
using ComPWA::Physics::SubSystem;

void EvtGenIF::addResonance(std::string name, double m0, double g0, double spin,
                            SubSystem subsys) {
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

  evtPars[std::string(name + "_mass")] = std::shared_ptr<ComPWA::FitParameter>(
      new ComPWA::FitParameter(name + "_mass", m0));
  evtPars.at(std::string(name + "_mass"))->fixParameter(false);
  evtPars[std::string(name + "_width")] = std::shared_ptr<ComPWA::FitParameter>(
      new ComPWA::FitParameter(name + "_width", g0));
  evtPars.at(std::string(name + "_width"))->fixParameter(true);
  // LOG(debug) << "EvtGenIF::addResonance finishes";
}

void EvtGenIF::addHeliResonance(boost::property_tree::ptree pt,
                                std::shared_ptr<PartList> partL) {
  // LOG(debug) << "EvtGenIF::addHeliResonance starts";

  SubSystem SubSys(pt);
  // std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD>
  // AngularDist;
  // std::shared_ptr<ComPWA::Physics::DecayDynamics::AbstractDynamicalFunction>
  //   DynamicFcn;
  std::pair<std::string, std::string> DecayProducts;
  std::pair<ComPWA::Spin, ComPWA::Spin> DecayHelicities;

  // int DataPosition;
  // DataPosition =
  //    3 * std::dynamic_pointer_cast<HelicityKinematics>(kin)->dataID(SubSys);

  std::string resoName = pt.get<std::string>("<xmlattr>.Name", "empty");
  double magnitude;
  double phase;
  for (const auto &v : pt.get_child("")) {
    if (v.first == "Parameter") {
      if (v.second.get<std::string>("<xmlattr>.Type") == "Magnitude") {
        magnitude = (std::make_shared<FitParameter>(v.second))->value();
      }
      if (v.second.get<std::string>("<xmlattr>.Type") == "Phase") {
        phase = (std::make_shared<FitParameter>(v.second))->value();
      }
    }
  }
  // LOG(debug) << "EvtGenIF::addHeliResonance decay";

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  auto partItr = partL->find(name);
  if (partItr == partL->end())
    throw std::runtime_error("EvtGenIF::addHeliResonance() | Particle " + name +
                             " not found in list!");
  ComPWA::Spin J = partItr->second.GetSpinQuantumNumber("Spin");
  ComPWA::Spin mu(pt.get<double>("DecayParticle.<xmlattr>.Helicity"));

  // LOG(debug) << "EvtGenIF::addHeliResonance decay products";
  // Read name and helicities from decay products
  auto decayProducts = pt.get_child("DecayProducts");
  if (decayProducts.size() != 2)
    throw boost::property_tree::ptree_error(
        "EvtGenIF::addHeliResonance() | Expect exactly two decay products (" +
        std::to_string(decayProducts.size()) + " given)!");

  auto p = decayProducts.begin();
  DecayProducts.first = p->second.get<std::string>("<xmlattr>.Name");
  DecayHelicities.first =
      ComPWA::Spin(p->second.get<int>("<xmlattr>.Helicity"));
  ++p;
  DecayProducts.second = p->second.get<std::string>("<xmlattr>.Name");
  DecayHelicities.second =
      ComPWA::Spin(p->second.get<double>("<xmlattr>.Helicity"));

  // LOG(debug) << "EvtGenIF::addHeliResonance two-body decay";
  // Two-body decay
  if (SubSys.getFinalStates().size() == 2) {
    // Create WignerD object
    // AngularDist = std::make_shared<HelicityFormalism::AmpWignerD>(
    //     J, mu, DecayHelicities.first - DecayHelicities.second);

    auto partProp = partL->find(name)->second;
    std::string decayType = partProp.GetDecayType();

    if (decayType == "stable") {
      throw std::runtime_error(
          "EvtGenIF::addHeliResonance() | Stable particle is "
          "given as mother particle of a decay. Makes no "
          "sense!");
    } else if (decayType == "relativisticBreitWigner") {
      LOG(DEBUG) << "EvtGenIF::addHeliResonance add resonance";
      double width = 0.1;
      for (const auto &v : partProp.GetDecayInfo().get_child("")) {
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
      addResonance(name, partItr->second.GetMass(), width, J, SubSys);
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

void EvtGenIF::addResonances(boost::property_tree::ptree pt,
                             std::shared_ptr<DalitzKinematics> kin,
                             std::shared_ptr<PartList> partL) {
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
              if (x.second.get<std::string>("<xmlattr>.Type") == "Magnitude")
                double magnitude =
                    (std::make_shared<FitParameter>(x.second))->value();
              if (x.second.get<std::string>("<xmlattr>.Type") == "Phase")
                double phase =
                    (std::make_shared<FitParameter>(x.second))->value();
            } else if (x.first == "PartialAmplitude" &&
                       x.second.get<std::string>("<xmlattr>.Class") ==
                           "HelicityDecay") {
              LOG(DEBUG) << "EvtGenIF::addResoances add Heli-Resonance";
              addHeliResonance(x.second, partL);
              //} else if (x.first == "PartialAmplitude" &&
              //           x.second.get<std::string>("<xmlattr>.Class") ==
              //               "NonResonant") {
              //  addPartialAmplitude(
              //      std::make_shared<ComPWA::Physics::NonResonant>(partL, kin,
              //      x.second));
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

double EvtGenIF::intensity(const ComPWA::DataPoint &point) const {

  // We have to get around the constness of the interface definition.
  // std::vector<std::vector<double>> parameters(Parameters);

  std::vector<double> normValues(NormalizationValues);

  double result = 0;
  // std::cout << "EvtGenDalitzResos: " << Resos.size() << std::endl;
  for (unsigned int i = 0; i < Resos.size(); ++i) {
    EvtDalitzPoint pnt(EvtGenIF::transformToEvt(point));
    EvtDalitzReso tmp = Resos.at(i);
    std::string name = tmp.get_Name();
    // if(i+1<evtPars.size()){
    tmp.set_Mass(evtPars.at(std::string(name + "_mass"))->value());
    tmp.set_Gamma(evtPars.at(std::string(name + "_width"))->value());
    //}
    result += abs2(tmp.evaluate(pnt)); // * normValues.at(i);
  }
  // std::cout << "Result: " << result << std::endl;

  // const_cast<std::vector<std::vector<double>> &>(Parameters) = parameters;
  const_cast<std::vector<double> &>(NormalizationValues) = normValues;

  assert(!std::isnan(result) &&
         "IncoherentIntensity::Intensity() | Result is NaN!");
  assert(!std::isinf(result) &&
         "IncoherentIntensity::Intensity() | Result is inf!");

  return (strength() * result);
}

void EvtGenIF::parameters(ComPWA::ParameterList &list) {
  // Strength = list.addUniqueParameter(Strength);

  for (auto i : evtPars) {
    list.addParameter(i.second);
  }
}

void EvtGenIF::updateParameters(const ParameterList &list) {
  std::shared_ptr<FitParameter> p;
  try {
    p = FindParameter(Strength->name(), list);
  } catch (std::exception &ex) {
  }
  if (p)
    Strength->updateParameter(p);
  // for (auto i : Intensities)
  //  i->updateParameters(list);

  for (auto i : evtPars) {
    std::string name = i.first;
    std::shared_ptr<FitParameter> tmp;
    tmp = FindParameter(name, list);
    if (tmp) {
      evtPars[name]->setValue(tmp->value());
    }
  }

  return;
}

boost::property_tree::ptree EvtGenIF::save() const {
  boost::property_tree::ptree pt;
  // pt.put<std::string>("<xmlattr>.Name", name());
  // pt.add_child("Parameter", Strength->save());
  // pt.put("Parameter.<xmlattr>.Type", "Strength");
  // for (auto i : Intensities)
  //  pt.add_child("CoherentIntensity", i->save());

  return pt;
}

std::shared_ptr<ComPWA::AmpIntensity> EvtGenIF::component(std::string name) {

  return nullptr;
}

std::shared_ptr<ComPWA::AmpIntensity> EvtGenIF::component(std::string name,
    std::string resName, std::string daug1Name, std::string daug2Name,
    int L, int S) {

  return nullptr;
}

std::shared_ptr<ComPWA::FunctionTree>
EvtGenIF::tree(std::shared_ptr<Kinematics> kin,
               const ComPWA::ParameterList &sample,
               const ComPWA::ParameterList &phspSample,
               const ComPWA::ParameterList &toySample, unsigned int nEvtVar,
               std::string suffix) {

  return nullptr;
}
