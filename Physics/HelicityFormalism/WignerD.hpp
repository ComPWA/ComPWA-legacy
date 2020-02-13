// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_HELICITY_FORMALISM_WIGNERD_HPP_
#define COMPWA_PHYSICS_HELICITY_FORMALISM_WIGNERD_HPP_

#include "Core/FunctionTree/Functions.hpp"
#include "Core/FunctionTree/TreeNode.hpp"

#include "ThirdParty/qft++/include/qft++/WignerD.h"

namespace ComPWA {

struct DataPoint;

namespace Physics {
namespace HelicityFormalism {

///
/// Angular distribution based on WignerD functions
///
namespace WignerD {

inline double dynamicalFunction(double J, double muPrime, double mu,
                                double beta) {

  if ((double)J == 0)
    return 1.0;

  assert(!std::isnan(beta));
  assert(std::cos(beta) <= 1 && std::cos(beta) >= -1);

  double result = QFT::Wigner_d(J, muPrime, mu, beta);
  assert(!std::isnan(result));

  double pi4 = M_PI * 4.0;
  double norm = std::sqrt((2.0 * J + 1) / pi4);

  return norm * result;
}

inline std::complex<double> dynamicalFunction(double J, double muPrime,
                                              double mu, double alpha,
                                              double beta, double gamma) {
  if ((double)J == 0)
    return std::complex<double>(1.0, 0);

  assert(!std::isnan(alpha));
  assert(!std::isnan(beta));
  assert(!std::isnan(gamma));

  std::complex<double> i(0, 1);

  double tmp = WignerD::dynamicalFunction(J, muPrime, mu, beta);
  std::complex<double> result =
      tmp * std::exp(-i * (muPrime * alpha + mu * gamma));

  assert(!std::isnan(result.real()));
  assert(!std::isnan(result.imag()));

  return result;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode> createFunctionTree(
    double J, double MuPrime, double Mu,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>> Theta,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>> Phi);
} // namespace WignerD

class WignerDStrategy : public ComPWA::FunctionTree::Strategy {
public:
  WignerDStrategy(std::string Name = "")
      : Strategy(ComPWA::FunctionTree::ParType::MCOMPLEX, "WignerD" + Name) {}

  virtual void execute(ComPWA::FunctionTree::ParameterList &paras,
                       std::shared_ptr<ComPWA::FunctionTree::Parameter> &out);
};

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA

#endif
