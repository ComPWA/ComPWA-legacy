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
  return IntensityNoEff(point) * _eff->Evaluate(point);
}

double CoherentIntensity::IntensityNoEff(const dataPoint &point) const {
  std::complex<double> result(0., 0.);
  for (auto i : _seqDecays)
    result += i->Evaluate(point);
  return GetStrengthValue() * std::norm(result) * GetNormalization();
};

std::shared_ptr<CoherentIntensity>
CoherentIntensity::Factory(const boost::property_tree::ptree &pt) {
  LOG(trace) << " CoherentIntensity::Factory() | Construction....";
  auto obj = std::make_shared<CoherentIntensity>();
  obj->SetName(pt.get<std::string>("<xmlattr>.Name"));

  //  boost::property_tree::xml_writer_settings<char> settings('\t', 1);
  //  write_xml(std::cout,pt);

  auto ptCh = pt.get_child_optional("Strength");
  if (ptCh) {
    auto strength = ComPWA::DoubleParameterFactory(ptCh.get());
    obj->SetStrength(std::make_shared<DoubleParameter>(strength));
  } else {
    obj->SetStrength(std::make_shared<ComPWA::DoubleParameter>("", 1.0));
  }

  for (const auto &v : pt.get_child("")) {
    if (v.first == "Amplitude")
      obj->Add(
          ComPWA::Physics::HelicityFormalism::SequentialTwoBodyDecay::Factory(
              v.second));
  }
  return obj;
}

boost::property_tree::ptree
CoherentIntensity::Save(std::shared_ptr<CoherentIntensity> obj) {

  boost::property_tree::ptree pt;
  pt.put<std::string>("<xmlattr>.Name", obj->GetName());
  pt.add_child("Strength",
               ComPWA::DoubleParameterSave(*obj->GetStrength().get()));
  for (auto i : obj->GetDecays()) {
    pt.add_child("Amplitude", SequentialTwoBodyDecay::Save(i));
  }
  return pt;
}

//! Getter function for basic amp tree
std::shared_ptr<ComPWA::FunctionTree> CoherentIntensity::GetTree(
    ComPWA::ParameterList &sample, ComPWA::ParameterList &phspSample,
    ComPWA::ParameterList &toySample, std::string suffix) {

  unsigned int effId = Kinematics::Instance()->GetNVars();
  unsigned int weightId = Kinematics::Instance()->GetNVars() + 1;
  int phspSampleSize = phspSample.GetMultiDouble(0)->GetNValues();

  std::shared_ptr<MultiDouble> weightPhsp = phspSample.GetMultiDouble(weightId);
  double sumWeights =
      std::accumulate(weightPhsp->Begin(), weightPhsp->End(), 0.0);
  std::shared_ptr<MultiDouble> eff = phspSample.GetMultiDouble(effId);

  std::shared_ptr<FunctionTree> tr(new FunctionTree());

  tr->createHead("CoherentIntens(" + GetName() + ")" + suffix,
                 std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  tr->createLeaf("Strength", _strength,
                 "CoherentIntens(" + GetName() + ")" + suffix);
  tr->createNode("SumSquared",
                 std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                 "CoherentIntens(" + GetName() + ")" + suffix);
  tr->insertTree(setupBasicTree(sample, toySample), "SumSquared");

  // Normalization
  // create a new FunctionTree to make sure that nodes with the same name do
  // not interfere
  std::shared_ptr<FunctionTree> trNorm(new FunctionTree());
  trNorm->createHead(
      "CoherentIntensNorm",
      std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE))); // 1/normLH
  trNorm->createNode("NormFactor",
                     std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
                     "CoherentIntensNorm");
  trNorm->createLeaf("PhspVolume", Kinematics::Instance()->GetPhspVolume(),
                     "NormFactor");
  trNorm->createLeaf("InvNmc", 1 / ((double)sumWeights), "NormFactor");
  trNorm->createNode("SumAmp",
                     std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                     "NormFactor");
  trNorm->createNode("IntensPhspEff",
                     std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)),
                     "SumAmp", phspSampleSize, false);
  trNorm->createLeaf("Eff", eff, "IntensPhspEff");
  trNorm->createLeaf("WeightPhsp", weightPhsp, "IntensPhspEff");
  trNorm->createNode("IntensPhsp",
                     std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                     "IntensPhspEff", phspSampleSize,
                     false); //|T_{ev}|^2
  trNorm->insertTree(setupBasicTree(phspSample, toySample, "_norm"),
                     "IntensPhsp");

  tr->insertTree(trNorm, "CoherentIntens(" + GetName() + ")" + suffix);
  return tr;
}

std::shared_ptr<FunctionTree>
CoherentIntensity::setupBasicTree(ParameterList &sample,
                                  ParameterList &phspSample,
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
  list.AddParameter(GetStrength());
  for (auto i : _seqDecays) {
    i->GetParameters(list);
  }
}

double CoherentIntensity::GetNormalization() const {
  if (_integral)
    return 1 / _integral;

  // Check if parameters were modified
  ParameterList list;
  const_cast<CoherentIntensity *>(this)->GetParameters(list);
  if (const_cast<ParameterList &>(_currentParList) == list && _integral > 0 &&
      _maxIntens > 0)
    return 1 / _integral;

  const_cast<ParameterList &>(_currentParList).DeepCopy(list);
  // Parameters have changed. We have to recalculate the normalization.
  const_cast<double &>(_integral) = Integral();

  return 1 / _integral;
}

double CoherentIntensity::Integral() const {

  if (!_phspSample->size()) {
    LOG(debug)
        << "CoherentIntensity::Integral() | Integral can not be calculated "
           "since no phsp sample is set. Set a sample using "
           "SetPhspSamples(phspSample, toySample)!";
    return 1.0;
  }

  double sumIntens = 0;
  double maxVal = 0;
  for (auto i : *_phspSample.get()) {
    std::complex<double> result(0., 0.);
    for (auto d : _seqDecays)
      result += d->Evaluate(i);
    double intens = std::norm(result);

    if (intens > maxVal)
      maxVal = intens;

    sumIntens += intens;
  }

  double phspVol = Kinematics::Instance()->GetPhspVolume();
  double integral = (sumIntens * phspVol / _phspSample->size());
  LOG(trace) << "CoherentIntensity::Integral() | Integral is " << integral
             << " and the maximum value of intensity is " << maxVal << ".";

  const_cast<double &>(_maxIntens) = maxVal;
  return integral;
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
