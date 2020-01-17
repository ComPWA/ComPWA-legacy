// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "RelativisticBreitWigner.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

using ComPWA::FunctionTree::FitParameter;
using ComPWA::FunctionTree::Parameter;
using ComPWA::FunctionTree::ParameterList;
using ComPWA::FunctionTree::TreeNode;
using ComPWA::FunctionTree::Value;

std::shared_ptr<TreeNode> RelativisticBreitWigner::createFunctionTree(
    RelativisticBreitWigner::InputInfo Params,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquared) {
  size_t sampleSize = InvMassSquared->values().size();

  BreitWignerFunction BWFunction(
      RelativisticBreitWigner::relativisticBreitWigner);
  if (Params.Type == "relativisticBreitWigner") {
  } else if (Params.Type == "relativisticBreitWignerAC") {
    BWFunction = RelativisticBreitWigner::relativisticBreitWignerAnalyticCont;
  } else {
    LOG(INFO) << "Relativistic BreitWigner of type " << Params.Type
              << " is unknown. Using standard defintion as fallback!";
  }

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<TreeNode>(
      MComplex("", sampleSize), std::make_shared<BreitWignerStrategy>(
                                    Params.FormFactorFunctor, BWFunction));

  tr->addNodes({createLeaf(Params.Mass), createLeaf(Params.Width),
                createLeaf((int)Params.L, "L"), createLeaf(Params.MesonRadius),
                createLeaf(Params.DaughterMasses.first),
                createLeaf(Params.DaughterMasses.second),
                createLeaf(InvMassSquared)});

  return tr;
};

void BreitWignerStrategy::execute(ParameterList &paras,
                                  std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter(
        "BreitWignerStrat::execute() | Parameter type mismatch!");

#ifndef NDEBUG
  // Check parameter type
  if (checkType != out->type())
    throw(
        WrongParType("BreitWignerStrat::execute() | "
                     "Output parameter is of type " +
                     std::string(ComPWA::FunctionTree::ParNames[out->type()]) +
                     " and conflicts with expected type " +
                     std::string(ComPWA::FunctionTree::ParNames[checkType])));

  // How many parameters do we expect?
  size_t check_nInt = 1;
  size_t nInt = paras.intValues().size();
  size_t check_nDouble = 5;
  size_t nDouble = paras.doubleValues().size();
  nDouble += paras.doubleParameters().size();
  size_t check_nComplex = 0;
  size_t nComplex = paras.complexValues().size();
  size_t check_nMInteger = 0;
  size_t nMInteger = paras.mIntValues().size();
  size_t check_nMDouble = 1;
  size_t nMDouble = paras.mDoubleValues().size();
  size_t check_nMComplex = 0;
  size_t nMComplex = paras.mComplexValues().size();

  // Check size of parameter list
  if (nInt != check_nInt)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of IntParameters does not match: " +
                       std::to_string(nInt) + " given but " +
                       std::to_string(check_nInt) + " expected."));
  if (nDouble != check_nDouble)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of FitParameters does not match: " +
                       std::to_string(nDouble) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (nComplex != check_nComplex)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(nComplex) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (nMInteger != check_nMInteger)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of MultiInt does not match: " +
                       std::to_string(nMInteger) + " given but " +
                       std::to_string(check_nMInteger) + " expected."));
  if (nMDouble != check_nMDouble)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of MultiDoubles does not match: " +
                       std::to_string(nMDouble) + " given but " +
                       std::to_string(check_nMDouble) + " expected."));
  if (nMComplex != check_nMComplex)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of MultiComplexes does not match: " +
                       std::to_string(nMComplex) + " given but " +
                       std::to_string(check_nMComplex) + " expected."));
#endif

  size_t n = paras.mDoubleValue(0)->values().size();
  if (!out)
    out = ComPWA::FunctionTree::MComplex("", n);
  auto par =
      std::static_pointer_cast<Value<std::vector<std::complex<double>>>>(out);
  auto &results = par->values(); // reference
  if (results.size() != n) {
    results.resize(n);
  }
  // Get parameters from ParameterList:
  // We use the same order of the parameters as was used during tree
  // construction.
  double m0 = paras.doubleParameter(0)->value();
  double Gamma0 = paras.doubleParameter(1)->value();
  double MesonRadius = paras.doubleParameter(2)->value();
  unsigned int orbitL = paras.intValue(0)->value();
  double ma = paras.doubleParameter(3)->value();
  double mb = paras.doubleParameter(4)->value();

  std::transform(paras.mDoubleValue(0)->values().begin(),
                 paras.mDoubleValue(0)->values().end(), results.begin(),
                 [&](double s) {
                   return BWFunction(s, m0, ma, mb, Gamma0, orbitL, MesonRadius,
                                     FormFactorFunctor);
                 });
}

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA
