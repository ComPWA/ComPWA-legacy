//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//--------------------------------------------------------------------------------
#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/Efficiency.hpp"
#include "Physics/HelicityFormalism/CoherentIntensity.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

double CoherentIntensity::Intensity(const dataPoint &point) const {
  std::complex<double> result(0., 0.);
  for (auto i : _seqDecays)
    result += i->Evaluate(point);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return Strength() * std::norm(result) * _eff->Evaluate(point);
};

std::shared_ptr<CoherentIntensity>
CoherentIntensity::Factory(const boost::property_tree::ptree &pt) {
  LOG(trace) << " CoherentIntensity::Factory() | Construction....";
  auto obj = std::make_shared<CoherentIntensity>();
  obj->_name=(pt.get<std::string>("<xmlattr>.Name"));

  //  boost::property_tree::xml_writer_settings<char> settings('\t', 1);
  //  write_xml(std::cout,pt);

  auto ptCh = pt.get_child_optional("Strength");
  if (ptCh) {
    auto strength = ComPWA::DoubleParameterFactory(ptCh.get());
    obj->_strength=(std::make_shared<DoubleParameter>(strength));
  } else {
    obj->_strength=(
        std::make_shared<ComPWA::DoubleParameter>("", 1.0));
  }

  for (const auto &v : pt.get_child("")) {
    if (v.first == "Amplitude")
      obj->AddAmplitude(
          ComPWA::Physics::HelicityFormalism::SequentialTwoBodyDecay::Factory(
              v.second));
  }
  return obj;
}

boost::property_tree::ptree
CoherentIntensity::Save(std::shared_ptr<CoherentIntensity> obj) {

  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->Name());
  pt.add_child("Strength",
               ComPWA::DoubleParameterSave(*obj->_strength.get()));
  for (auto i : obj->GetAmplitudes()) {
    pt.add_child("Amplitude", SequentialTwoBodyDecay::Save(i));
  }
  return pt;
}

std::shared_ptr<AmpIntensity>
CoherentIntensity::GetComponent(std::string name) {

  // The whole object?
  if (name == Name()) {
    LOG(error) << "CoherentIntensity::GetComponent() | You're requesting the "
                  "full object! So just copy it!";
    return std::shared_ptr<AmpIntensity>();
  }

  bool found = false;
  // Do we want to have a combination of coherentintensities?
  std::vector<std::string> splitNames = splitString(name);
  auto icIn = std::shared_ptr<AmpIntensity>(this->Clone(name));
  icIn->Reset(); //delete all existing amplitudes
  for (auto i : splitNames) {
    for (auto j : _seqDecays) {
      if (i == j->GetName()) {
        std::dynamic_pointer_cast<CoherentIntensity>(icIn)->AddAmplitude(j);
        found = true;
      }
    }
  }

  // Nothing found
  if (!found) {
    throw std::runtime_error(
        "CoherentIntensity::GetComponent() | Component " + name +
        " could not be found in CoherentIntensity " + Name() + ".");
  }
  
  return icIn;
}

//! Getter function for basic amp tree
std::shared_ptr<ComPWA::FunctionTree>
CoherentIntensity::GetTree(const ComPWA::ParameterList &sample,
                           const ComPWA::ParameterList &phspSample,
                           const ComPWA::ParameterList &toySample,
                           std::string suffix) {

  unsigned int effId = Kinematics::Instance()->GetNVars();
  unsigned int weightId = Kinematics::Instance()->GetNVars() + 1;
  int phspSampleSize = phspSample.GetMultiDouble(0)->GetNValues();

  std::shared_ptr<MultiDouble> weightPhsp = phspSample.GetMultiDouble(weightId);
  double sumWeights =
      std::accumulate(weightPhsp->Begin(), weightPhsp->End(), 0.0);
  std::shared_ptr<MultiDouble> eff = phspSample.GetMultiDouble(effId);

  std::shared_ptr<FunctionTree> tr(new FunctionTree());

  tr->createHead("CoherentIntensity(" + Name() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  tr->createLeaf("Strength", _strength,
                 "CoherentIntensity(" + Name() + ")" + suffix);
  tr->createNode("SumSquared",
                 std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                 "CoherentIntensity(" + Name() + ")" + suffix);
  tr->insertTree(setupBasicTree(sample, toySample), "SumSquared");

  // Normalization
  // create a new FunctionTree to make sure that nodes with the same name do
  // not interfere
  std::shared_ptr<FunctionTree> trNorm(new FunctionTree());
  trNorm->createHead("Normalization(" + Name() + ")" + suffix,
                     std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)));
  trNorm->createNode("Integral",
                     std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
                     "Normalization(" + Name() + ")" + suffix);
  trNorm->createLeaf("PhspVolume", Kinematics::Instance()->GetPhspVolume(),
                     "Integral");
  trNorm->createLeaf("InverseSampleWeights", 1 / ((double)sumWeights),
                     "Integral");
  trNorm->createNode("Sum",
                     std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                     "Integral");
  trNorm->createNode("IntensityWeighted",
                     std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)),
                     "Sum", phspSampleSize, false);
  trNorm->createLeaf("Efficiency", eff, "IntensityWeighted");
  trNorm->createLeaf("EventWeight", weightPhsp, "IntensityWeighted");
  trNorm->createNode("Intensity",
                     std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                     "IntensityWeighted", phspSampleSize,
                     false); //|T_{ev}|^2
  trNorm->insertTree(setupBasicTree(phspSample, toySample, "_norm"),
                     "Intensity");

  tr->insertTree(trNorm, "CoherentIntensity(" + Name() + ")" + suffix);
  return tr;
}

std::shared_ptr<FunctionTree>
CoherentIntensity::setupBasicTree(const ParameterList &sample,
                                  const ParameterList &phspSample,
                                  std::string suffix) const {

  int sampleSize = sample.GetMultiDouble(0)->GetNValues();
  int phspSampleSize = phspSample.GetMultiDouble(0)->GetNValues();

  if (sampleSize == 0) {
    LOG(error) << "AmpSumIntensity::setupBasicTree() | "
                  "Data sample empty!";
    return std::shared_ptr<FunctionTree>();
  }
  if (phspSampleSize == 0) {
    LOG(error) << "AmpSumIntensity::setupBasicTree() | "
                  "Phsp sample empty!";
    return std::shared_ptr<FunctionTree>();
  }

  //------------Setup Tree---------------------
  std::shared_ptr<FunctionTree> newTree(new FunctionTree());

  //----Strategies needed
  std::shared_ptr<AddAll> maddStrat(new AddAll(ParType::MCOMPLEX));

  newTree->createHead("SumOfAmplitudes" + suffix, maddStrat, sampleSize);

  for (auto i : _seqDecays) {
    std::shared_ptr<FunctionTree> resTree = i->GetTree(sample, phspSample, "");
    if (!resTree->sanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Resonance tree didn't pass sanity check!");
    resTree->recalculate();
    newTree->insertTree(resTree, "SumOfAmplitudes" + suffix);
  }

  LOG(debug) << "AmpSumIntensity::setupBasicTree(): tree constructed!!";
  return newTree;
}

void CoherentIntensity::GetParameters(ComPWA::ParameterList &list) {
  list.AddParameter(_strength);
  for (auto i : _seqDecays) {
    i->GetParameters(list);
  }
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
