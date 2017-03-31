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
  return GetStrengthValue() * std::norm(result);
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
  tr->createNode("N", std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)),
                 "CoherentIntens(" + GetName() + ")" + suffix); // 1/normLH
  // normLH = phspVolume/N_{mc} |T_{evPHSP}|^2
  tr->createNode("normFactor",
                 std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)), "N");
  // sumAmp = \sum_{evPHSP} |T_{evPHSP}|^2
  tr->createNode("sumAmp",
                 std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                 "normFactor");
  tr->createLeaf("phspVolume", Kinematics::Instance()->GetPhspVolume(),
                 "normFactor");
  tr->createLeaf("InvNmc", 1 / ((double)sumWeights), "normFactor");
  tr->createNode("IntensPhspEff",
                 std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)),
                 "sumAmp", phspSampleSize,
                 false);                       //|T_{ev}|^2
  tr->createLeaf("eff", eff, "IntensPhspEff"); // efficiency
  tr->createLeaf("weightPhsp", weightPhsp, "IntensPhspEff");
  tr->createNode("IntensPhsp",
                 std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                 "IntensPhspEff", phspSampleSize,
                 false); //|T_{ev}|^2
  tr->insertTree(setupBasicTree(phspSample, toySample, "_norm"), "IntensPhsp");

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
    std::shared_ptr<FunctionTree> resTree =
        i->GetTree(sample, phspSample, phspSample, "");
    if (!resTree->sanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Resonance tree didn't pass sanity check!");
    resTree->recalculate();
    newTree->insertTree(resTree, "SumOfAmplitudes" + suffix);
  }

  LOG(debug) << "AmpSumIntensity::setupBasicTree(): tree constructed!!";
  return newTree;
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
