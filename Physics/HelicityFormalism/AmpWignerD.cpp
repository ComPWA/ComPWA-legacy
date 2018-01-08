// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
#include "Physics/HelicityFormalism/AmpWignerD.hpp"
#include "Physics/qft++/WignerD.h"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

AmpWignerD::AmpWignerD(ComPWA::Spin spin, ComPWA::Spin mu, ComPWA::Spin muPrime)
    : J(spin), Mu(mu), MuPrime(muPrime) {}

double AmpWignerD::evaluate(const DataPoint &point, int pos1, int pos2) const {
  if ((double)J == 0)
    return 1.0;
  return dynamicalFunction(J, Mu, MuPrime, point.value(pos1));
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
  //  double norm = 1/sqrt(2*J.GetSpin()+1);
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

std::shared_ptr<FunctionTree> AmpWignerD::tree(const ParameterList &sample,
                                                  int posTheta, int posPhi,
                                                  std::string suffix) {

  size_t n = sample.mDoubleValue(0)->values().size();

  // in case of spin zero do not explicitly include the WignerD
  if ((double)J == 0)
    return std::make_shared<FunctionTree>("WignerD" + suffix,
                                          MInteger("", n, 1));

  auto tr = std::make_shared<FunctionTree>(
      "WignerD" + suffix, MDouble("", n),
      std::make_shared<WignerDStrategy>("WignerD" + suffix));

  tr->createLeaf("spin", (double)J, "WignerD" + suffix); // spin
  tr->createLeaf("m", (double)Mu, "WignerD" + suffix);      // OutSpin 1
  tr->createLeaf("n", (double)(MuPrime), "WignerD" + suffix); // OutSpin 2
  tr->createLeaf("data_cosTheta[" + std::to_string(posTheta) + "]",
                 sample.mDoubleValue(posTheta), "WignerD" + suffix);
  tr->createLeaf("data_phi[" + std::to_string(posPhi) + "]",
                 sample.mDoubleValue(posPhi), "WignerD" + suffix);

  return tr;
}

void WignerDStrategy::execute(ParameterList &paras,
                              std::shared_ptr<Parameter> &out) {
#ifndef NDEBUG
  if (out && checkType != out->type()) {
    throw(WrongParType(std::string("Output Type ") + ParNames[out->type()] +
                       std::string(" conflicts expected type ") +
                       ParNames[checkType] + std::string(" of ") + name +
                       " Wigner strat"));
  }
#endif

  double J = paras.doubleValue(0)->value();
  double mu = paras.doubleValue(1)->value();
  double muPrime = paras.doubleValue(2)->value();

  auto data = paras.mDoubleValue(0);

  size_t n = data->values().size();
  if (!out)
    out = MDouble("", n);
  auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
  auto &results = par->values(); // reference

  for (unsigned int ele = 0; ele < n; ele++) {
    try {
      results.at(ele) = AmpWignerD::dynamicalFunction(
          J, mu, muPrime, data->values().at(ele));
    } catch (std::exception &ex) {
      LOG(error) << "WignerDStrategy::execute() | " << ex.what();
      throw std::runtime_error("WignerDStrategy::execute() | "
                               "Evaluation of dynamical function failed!");
    }
  } // end element loop
}

} // ns::HelicityFormalism
} // ns::Physics
} // ns::ComPWA
