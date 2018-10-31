// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/Efficiency.hpp"
#include "Physics/CoherentIntensity.hpp"
#include "Tools/Integration.hpp"

using namespace ComPWA;
using namespace ComPWA::Physics;
using namespace ComPWA::Physics::HelicityFormalism;

CoherentIntensity::CoherentIntensity(std::shared_ptr<ComPWA::PartList> partL,
                  std::shared_ptr<ComPWA::Kinematics> kin,
                  const boost::property_tree::ptree &pt) : PhspVolume(1.0) {
  load(partL, kin, pt);
}

double CoherentIntensity::intensity(const DataPoint &point) const {
  std::complex<double> result(0., 0.);
  for (auto i : Amplitudes)
    result += i->evaluate(point);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return strength() * std::norm(result) * Eff->evaluate(point);
};

void CoherentIntensity::load(std::shared_ptr<PartList> partL,
                             std::shared_ptr<Kinematics> kin,
                             const boost::property_tree::ptree &pt) {
  LOG(TRACE) << "CoherentIntensity::load() | Construction....";
  if (pt.get<std::string>("<xmlattr>.Class") != "Coherent")
    throw BadConfig("CoherentIntensity::Factory() | Property tree seems to "
                    "not containt a configuration for an "
                    "CoherentIntensity!");
  Name = (pt.get<std::string>("<xmlattr>.Name"));
  Strength = std::make_shared<ComPWA::FitParameter>("Strength_"+Name, 1.0);
  setPhspVolume(kin->phspVolume());
  
  for (const auto &v : pt.get_child("")) {
    // Strength parameter
    if (v.first == "Parameter" &&
        v.second.get<std::string>("<xmlattr>.Type") == "Strength") {
      Strength = std::make_shared<FitParameter>(v.second);
    } else if (v.first == "Amplitude" &&
               v.second.get<std::string>("<xmlattr>.Class") ==
                   "SequentialPartialAmplitude") {
      addAmplitude(
          std::make_shared<SequentialPartialAmplitude>(partL, kin, v.second));
    } else {
      // ignored further settings. Should we throw an error?
    }
  }
}

boost::property_tree::ptree CoherentIntensity::save() const {
  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", name());
  pt.add_child("Parameter", Strength->save());
  pt.put("Parameter.<xmlattr>.Type", "Strength");

  for (auto i : Amplitudes)
    pt.add_child("Amplitude", i->save());
  
  return pt;
}

std::shared_ptr<ComPWA::AmpIntensity>
CoherentIntensity::component(std::string name) {

  // The whole object?
  if (name == Name) {
    //    LOG(ERROR) << "CoherentIntensity::component() | You're requesting the
    //    "
    //                  "full object! So just copy it!";
    //    return std::shared_ptr<AmpIntensity>();
    return shared_from_this();
  }

  bool found = false;
  // Do we want to have a combination of coherentintensities?
  std::vector<std::string> splitNames = splitString(name);
  auto icIn = std::make_shared<CoherentIntensity>(*this);
  icIn->setName(name);
  icIn->reset(); // delete all existing amplitudes
  for (auto i : splitNames) {
    for (auto j : Amplitudes) {
      if (i == j->name()) {
        icIn->addAmplitude(j);
        found = true;
      }
    }
  }

  // Nothing found
  if (!found) {
    throw std::runtime_error(
        "CoherentIntensity::component() | Component " + name +
        " could not be found in CoherentIntensity " + Name + ".");
  }

  return icIn;
}

std::shared_ptr<ComPWA::AmpIntensity>
CoherentIntensity::component(std::string name, std::string resName, 
    std::string daug1Name, std::string daug2Name, int L, int S) {
  std::string lStr = "L_" + std::to_string(L) + ".0";
  if (L < 0) lStr = "";
  std::string sStr = "S_" + std::to_string(S) + ".0";
  if (S < 0) sStr = "";
  std::string lsStr = lStr + "_" + sStr;
  if (name == "") {
    name = resName + "_to_" + daug1Name + "+" + daug2Name;
    if (L >= 0) name += "_" + lStr;
    if (S >= 0) name += "_" + sStr;
  }


  auto icIn = std::make_shared<CoherentIntensity>(*this);
  icIn->setName(name);
  icIn->reset(); // delete all existing amplitudes
  bool found(false);
  for (auto seqAmp : Amplitudes) {
    std::string seqAmpName = seqAmp->name();
    bool foundComponent(false);
    std::size_t begPos = 0;
    std::size_t endPos = 0;
    while (begPos != std::string::npos) {
      endPos = seqAmpName.find(";", begPos);
      //endPos == std::string::npos in substr means until end of the string
      std::string subName = seqAmpName.substr(begPos, endPos);
      begPos = endPos;
      if (std::string::npos == subName.find(resName))
        continue;
      if (daug1Name != "" && std::string::npos == subName.find(daug1Name))
        continue;
      if (daug2Name != "" && std::string::npos == subName.find(daug2Name))
        continue;
      if (L >= 0 && std::string::npos == subName.find(lStr))
        continue; 
      if (S >= 0 && std::string::npos == subName.find(sStr))
        continue;
      foundComponent = true;
      break;
    }
    if (!foundComponent)
      continue;
    found = true;
    icIn->addAmplitude(seqAmp);
  }

  // Nothing found, this could be true 
  // and it is not a error when this is called
  // in the IncoherentIntnesity::component()
  if (!found) {
    throw std::runtime_error(
        "CoherentIntensity::component() | Component " + name +
        " could not be found in CoherentIntensity " + Name + ".");
  }

  return icIn;
}

std::shared_ptr<ComPWA::FunctionTree>
CoherentIntensity::tree(std::shared_ptr<Kinematics> kin,
                        const ComPWA::ParameterList &sample,
                        const ComPWA::ParameterList &phspSample,
                        const ComPWA::ParameterList &toySample,
                        unsigned int nEvtVar, std::string suffix) {

  size_t n = sample.mDoubleValue(0)->values().size();

  auto tr = std::make_shared<FunctionTree>(
      "CoherentIntensity(" + name() + ")" + suffix, MDouble("", n),
      std::make_shared<MultAll>(ParType::MDOUBLE));

  tr->createLeaf("Strength", Strength,
                 "CoherentIntensity(" + name() + ")" + suffix);
  tr->createNode("SumSquared",
                 std::make_shared<AbsSquare>(ParType::MDOUBLE),
                 "CoherentIntensity(" + name() + ")" + suffix);
  tr->createNode("SumOfAmplitudes", MComplex("", n),
                 std::make_shared<AddAll>(ParType::MCOMPLEX),
                 "SumSquared");

  for (auto i : Amplitudes) {
    std::shared_ptr<ComPWA::FunctionTree> resTree =
        i->tree(kin, sample, toySample, "");
    if (!resTree->sanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Resonance tree didn't pass sanity check!");
    resTree->parameter();
    tr->insertTree(resTree, "SumOfAmplitudes" + suffix);
  }

  return tr;
}

void CoherentIntensity::parameters(ComPWA::ParameterList &list) {
  Strength = list.addUniqueParameter(Strength);
  for (auto i : Amplitudes) {
    i->parameters(list);
  }
}

void CoherentIntensity::updateParameters(const ParameterList &list) {
  std::shared_ptr<FitParameter> p;
  try {
    p = FindParameter(Strength->name(), list);
  } catch (std::exception &ex) {
  }
  if (p)
    Strength->updateParameter(p);

  for (auto i : Amplitudes)
    i->updateParameters(list);

  return;
}
