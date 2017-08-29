// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
#include "Physics/HelicityFormalism/AmpWignerD.hpp"
#include "Physics/qft++/WignerD.h"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

AmpWignerD::AmpWignerD(ComPWA::Spin spin, unsigned int mu, unsigned int muPrime)
    : _spin(spin), _mu(mu), _helicities(muPrime, 0) {}

double AmpWignerD::Evaluate(const dataPoint &point, int pos1, int pos2) const {
  if ((double)_spin == 0)
    return 1.0;
  return dynamicalFunction(_spin, _mu, (_helicities.first - _helicities.second),
                           point.GetValue(pos1));
}

double AmpWignerD::dynamicalFunction(ComPWA::Spin J, ComPWA::Spin mu,
                                     ComPWA::Spin muPrime, double cosTheta) {

  if ((double)J == 0)
    return 1.0;

  assert(!std::isnan(cosTheta));

  if (cosTheta > 1 || cosTheta < -1)
    throw std::runtime_error(
        "AmpWignerD::dynamicalFunction() | "
        "scattering angle out of range! Datapoint beyond phsp?");

  double result = QFT::Wigner_d(J, mu, muPrime, acos(cosTheta));
  assert(!std::isnan(result));

  // Not quite sure what the correct prefactor is in this case.
  //	double norm = 1/sqrt(2*J.GetSpin()+1);
  double norm = (2 * J.GetSpin() + 1);

  return norm * result;
}

std::complex<double>
AmpWignerD::dynamicalFunction(double cosAlpha, double cosBeta, double cosGamma,
                              ComPWA::Spin J, ComPWA::Spin mu,
                              ComPWA::Spin muPrime) {
  if ((double)J == 0)
    return 1.0;

  assert(!std::isnan(cosAlpha));
  assert(!std::isnan(cosBeta));
  assert(!std::isnan(cosGamma));

  std::complex<double> i(0, 1);

  double tmp = AmpWignerD::dynamicalFunction(J, mu, muPrime, cosBeta);
  std::complex<double> result =
      tmp * std::exp(-i * (mu.GetSpin() * acos(cosAlpha) +
                           muPrime.GetSpin() * acos(cosGamma)));

  assert(!std::isnan(result.real()));
  assert(!std::isnan(result.imag()));

  return result;
}

std::shared_ptr<AmpWignerD>
AmpWignerD::Factory(std::shared_ptr<PartList> partL,
                    const boost::property_tree::ptree &pt) {
  LOG(trace) << "AmpWignerD::Factory() | Construction....";
  auto obj = std::make_shared<AmpWignerD>();

  auto decayParticle = pt.get_child("DecayParticle");

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  ComPWA::Spin J = partL->find(name)->second.GetSpinQuantumNumber("Spin");
  obj->SetSpin(J);
  ComPWA::Spin mu(pt.get<double>("DecayParticle.<xmlattr>.Helicity"));
  obj->SetMu(mu);

  auto decayProducts = pt.get_child("SubSystem.DecayProducts");
  std::vector<ComPWA::Spin> vHelicity;
  for (auto i : decayProducts) {
    vHelicity.push_back(
        ComPWA::Spin(i.second.get<double>("<xmlattr>.Helicity")));
  }

  if (vHelicity.size() != 2)
    throw boost::property_tree::ptree_error(
        "AmpWignerD::Factory() | Expect exactly two decay products (" +
        std::to_string(decayProducts.size()) + " given)!");

  obj->SetMuPrime(vHelicity.at(0) - vHelicity.at(1));

  return obj;
}

std::shared_ptr<FunctionTree> AmpWignerD::GetTree(const ParameterList &sample,
                                                  int posTheta, int posPhi,
                                                  std::string suffix) {
  std::shared_ptr<FunctionTree> newTree(new FunctionTree());

  int sampleSize = sample.GetMultiDouble(0)->GetNValues();
  if ((double)_spin ==
      0) { // in case of spin zero do not explicitly include the WignerD
    std::shared_ptr<MultiUnsignedInteger> one(
        new MultiUnsignedInteger("", std::vector<unsigned int>(sampleSize, 1)));
    newTree->createLeaf("WignerD" + suffix, one, ""); // spin
    return newTree;
  }
  //----Strategies needed
  std::shared_ptr<WignerDStrategy> angdStrat(
      new WignerDStrategy("AngD" + suffix));
  newTree->createHead(
      "WignerD" + suffix,
      std::shared_ptr<WignerDStrategy>(new WignerDStrategy("WignerD" + suffix)),
      sampleSize);

  newTree->createLeaf("spin", (double)_spin, "WignerD" + suffix); // spin
  newTree->createLeaf("m", (double)_mu, "WignerD" + suffix);      // OutSpin 1
  newTree->createLeaf("n", (double)(_helicities.first - _helicities.second),
                      "WignerD" + suffix); // OutSpin 2
  newTree->createLeaf("data_cosTheta[" + std::to_string(posTheta) + "]",
                      sample.GetMultiDouble(posTheta), "WignerD" + suffix);
  newTree->createLeaf("data_phi[" + std::to_string(posPhi) + "]",
                      sample.GetMultiDouble(posPhi), "WignerD" + suffix);

  return newTree;
}

bool WignerDStrategy::execute(ParameterList &paras,
                              std::shared_ptr<AbsParameter> &out) {
#ifndef NDEBUG
  if (checkType != out->type()) {
    throw(WrongParType(std::string("Output Type ") + ParNames[out->type()] +
                       std::string(" conflicts expected type ") +
                       ParNames[checkType] + std::string(" of ") + name +
                       " Wigner strat"));
    return false;
  }
#endif

  double _inSpin = paras.GetDoubleParameter(0)->GetValue();
  double _outSpin1 = paras.GetDoubleParameter(1)->GetValue();
  double _outSpin2 = paras.GetDoubleParameter(2)->GetValue();

  std::shared_ptr<MultiDouble> _angle = paras.GetMultiDouble(0);

  std::vector<double> results(_angle->GetNValues(), 0.);
  for (unsigned int ele = 0; ele < _angle->GetNValues(); ele++) {
    try {
      results.at(ele) = AmpWignerD::dynamicalFunction(
          _inSpin, _outSpin1, _outSpin2, _angle->GetValue(ele));
    } catch (std::exception &ex) {
      LOG(error) << "WignerDStrategy::execute() | " << ex.what();
      throw std::runtime_error("WignerDStrategy::execute() | "
                               "Evaluation of dynamical function failed!");
    }
  } // end element loop
  out = std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(), results));

  return true;
}

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */
