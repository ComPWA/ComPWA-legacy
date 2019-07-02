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

AmpWignerD::AmpWignerD(ComPWA::Spin spin, ComPWA::Spin muPrime, ComPWA::Spin mu)
    : J(spin), MuPrime(muPrime), Mu(mu) {}

std::complex<double> AmpWignerD::evaluate(const DataPoint &point, int pos1,
                                          int pos2) const {
  if ((double)J == 0)
    return 1.0;
  double theta(point.KinematicVariableList[pos1]);
  double phi(point.KinematicVariableList[pos2]);
  // evaluating the Wigner D functions with gamma = 0 is crucial!
  // The reason is that that the calculated phi values in the kinematics class
  // follow this condition (except the top node, the phi values are
  // differences of the production and decay plane)
  return dynamicalFunction(J, MuPrime, Mu, phi, theta, 0);
}

double AmpWignerD::dynamicalFunction(ComPWA::Spin J, ComPWA::Spin muPrime,
                                     ComPWA::Spin mu, double beta) {

  if ((double)J == 0)
    return 1.0;

  assert(!std::isnan(beta));
  assert(std::cos(beta) <= 1 && std::cos(beta) >= -1);

  double result =
      QFT::Wigner_d(J.GetSpin(), muPrime.GetSpin(), mu.GetSpin(), beta);
  assert(!std::isnan(result));

  double pi4 = M_PI * 4.0;
  double norm = std::sqrt((2 * J.GetSpin() + 1) / pi4);

  return norm * result;
}

std::complex<double> AmpWignerD::dynamicalFunction(ComPWA::Spin J,
                                                   ComPWA::Spin muPrime,
                                                   ComPWA::Spin mu,
                                                   double alpha, double beta,
                                                   double gamma) {
  if ((double)J == 0)
    return std::complex<double>(1.0, 0);

  assert(!std::isnan(alpha));
  assert(!std::isnan(beta));
  assert(!std::isnan(gamma));

  std::complex<double> i(0, 1);

  double tmp = AmpWignerD::dynamicalFunction(J, muPrime, mu, beta);
  std::complex<double> result =
      tmp * std::exp(-i * (muPrime.GetSpin() * alpha + mu.GetSpin() * gamma));

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
  double muPrime = paras.doubleValue(1)->value();
  double mu = paras.doubleValue(2)->value();

  auto thetas = paras.mDoubleValue(0);
  auto phis = paras.mDoubleValue(1);

  size_t n = thetas->values().size();
  if (!out)
    out = MComplex("", n);
  auto par =
      std::static_pointer_cast<Value<std::vector<std::complex<double>>>>(out);
  auto &results = par->values(); // reference
  if (results.size() != n) {
    results.resize(n);
  }
  for (unsigned int ele = 0; ele < n; ele++) {
    try {
      results[ele] = AmpWignerD::dynamicalFunction(
          J, muPrime, mu, phis->values()[ele], thetas->values()[ele], 0.0);
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
