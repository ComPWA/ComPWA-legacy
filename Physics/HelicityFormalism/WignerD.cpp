// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "WignerD.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

using namespace ComPWA::FunctionTree;

std::shared_ptr<ComPWA::FunctionTree::TreeNode> WignerD::createFunctionTree(
    double J, double MuPrime, double Mu,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>> Theta,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>> Phi) {

  // in case of spin zero do not explicitly include the WignerD
  if ((double)J == 0)
    return FunctionTree::createLeaf(1, "WignerD");

  auto tr = std::make_shared<ComPWA::FunctionTree::TreeNode>(
      MComplex("", 0), std::make_shared<WignerDStrategy>());

  tr->addNodes({FunctionTree::createLeaf(J, "J"),
                FunctionTree::createLeaf(MuPrime, "mu'"),
                FunctionTree::createLeaf(Mu, "mu"),
                FunctionTree::createLeaf(Theta),
                FunctionTree::createLeaf(Phi)});

  return tr;
}

void WignerDStrategy::execute(
    ParameterList &paras,
    std::shared_ptr<ComPWA::FunctionTree::Parameter> &out) {
#ifndef NDEBUG
  if (out && checkType != out->type()) {
    throw(WrongParType(std::string("Output Type ") + ParNames[out->type()] +
                       std::string(" conflicts expected type ") +
                       ParNames[checkType] + std::string(" of ") + Name +
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
      results[ele] = WignerD::dynamicalFunction(
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
