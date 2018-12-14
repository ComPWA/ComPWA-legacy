// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>

#include "qft++/WignerD.h"

#include "Core/Event.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

AmpWignerD::AmpWignerD(ComPWA::Spin spin, ComPWA::Spin mu, ComPWA::Spin muPrime)
    : J(spin), Mu(mu), MuPrime(muPrime) {}

std::complex<double> AmpWignerD::evaluate(const DataPoint &point, int pos1,
                                          int pos2) const {
  if ((double)J == 0)
    return 1.0;
  double theta(point.KinematicVariableList[pos1]);
  double phi(point.KinematicVariableList[pos2]);
  return dynamicalFunction(J, Mu, MuPrime, theta, phi);
}

double AmpWignerD::dynamicalFunction(ComPWA::Spin J, ComPWA::Spin mu,
                                     ComPWA::Spin muPrime, double theta) {

  if ((double)J == 0)
    return 1.0;

  assert(!std::isnan(theta));
  assert(std::cos(theta) <= 1 && std::cos(theta) >= -1);

  double result =
      QFT::Wigner_d(J.GetSpin(), mu.GetSpin(), muPrime.GetSpin(), theta);
  assert(!std::isnan(result));

  // Not quite sure what the correct prefactor is in this case.
  //  double norm = 1/sqrt(2*J.GetSpin()+1);
  double norm = (2 * J.GetSpin() + 1);

  return norm * result;
}

std::complex<double> AmpWignerD::dynamicalFunction(ComPWA::Spin J,
                                                   ComPWA::Spin mu,
                                                   ComPWA::Spin muPrime,
                                                   double theta, double phi) {
  if ((double)J == 0)
    return 1.0;

  assert(!std::isnan(theta));
  assert(!std::isnan(phi));

  std::complex<double> i(0, 1);

  double tmp = AmpWignerD::dynamicalFunction(J, mu, muPrime, theta);
  std::complex<double> result =
      tmp * std::exp(i * (mu.GetSpin() - muPrime.GetSpin()) * phi);

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
      "WignerD" + suffix, MComplex("", n),
      std::make_shared<WignerDStrategy>("WignerD" + suffix));

  tr->createLeaf("spin", (double)J, "WignerD" + suffix);      // spin
  tr->createLeaf("m", (double)Mu, "WignerD" + suffix);        // OutSpin 1
  tr->createLeaf("n", (double)(MuPrime), "WignerD" + suffix); // OutSpin 2
  tr->createLeaf("data_theta[" + std::to_string(posTheta) + "]",
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

  auto thetas = paras.mDoubleValue(0);
  auto phis = paras.mDoubleValue(1);

  size_t n = thetas->values().size();
  if (!out)
    out = MComplex("", n);
  auto par =
      std::static_pointer_cast<Value<std::vector<std::complex<double>>>>(out);
  auto &results = par->values(); // reference

  for (unsigned int ele = 0; ele < n; ele++) {
    try {
      results[ele] = AmpWignerD::dynamicalFunction(
          J, mu, muPrime, thetas->values()[ele], phis->values()[ele]);
    } catch (std::exception &ex) {
      LOG(ERROR) << "WignerDStrategy::execute() | " << ex.what();
      throw std::runtime_error("WignerDStrategy::execute() | "
                               "Evaluation of dynamical function failed!");
    }
  } // end element loop
}

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA
