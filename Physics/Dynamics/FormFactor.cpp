// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "FormFactor.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

using ComPWA::FunctionTree::FitParameter;
using ComPWA::FunctionTree::Parameter;
using ComPWA::FunctionTree::ParameterList;
using ComPWA::FunctionTree::TreeNode;
using ComPWA::FunctionTree::Value;

std::shared_ptr<ComPWA::FunctionTree::TreeNode> createFunctionTree(
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquaredDaughter1,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquaredDaughter2,
    std::shared_ptr<ComPWA::FunctionTree::FitParameter> MesonRadius,
    unsigned int L, std::shared_ptr<FormFactor> FormFactorFunctor,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquared) {
  size_t sampleSize = InvMassSquared->values().size();

  auto ffTree = std::make_shared<TreeNode>(
      ComPWA::FunctionTree::MDouble("", sampleSize),
      std::make_shared<FormFactorStrategy>(FormFactorFunctor));
  // add L and FFType as double value leaf, since there is no int leaf
  ffTree->addNodes({FunctionTree::createLeaf((int)L, "L"),
                    FunctionTree::createLeaf(MesonRadius),
                    FunctionTree::createLeaf(InvMassSquared),
                    FunctionTree::createLeaf(InvMassSquaredDaughter1),
                    FunctionTree::createLeaf(InvMassSquaredDaughter2)});
  ffTree->parameter();

  return ffTree;
}

void FormFactorStrategy::execute(ParameterList &paras,
                                 std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("FormFactorStrat::execute() | Parameter type mismatch!");

#ifndef NDEBUG
  // Check parameter type
  if (checkType != out->type())
    throw(
        WrongParType("FormFactorStrat::execute() | "
                     "Output parameter is of type " +
                     std::string(ComPWA::FunctionTree::ParNames[out->type()]) +
                     " and conflicts with expected type " +
                     std::string(ComPWA::FunctionTree::ParNames[checkType])));

  // How many parameters do we expect?
  size_t check_nInt = 1;
  size_t nInt = paras.intValues().size();
  // MesonRadius, Daughter1Mass, Daughter2Mass
  size_t check_nDouble = 1;
  size_t nDouble = paras.doubleValues().size();
  nDouble += paras.doubleParameters().size();
  size_t check_nComplex = 0;
  size_t nComplex = paras.complexValues().size();
  size_t check_nMInteger = 0;
  size_t nMInteger = paras.mIntValues().size();
  // DataSample.mDoubleValue(pos) (mSq)
  size_t check_nMDouble = 3;
  size_t nMDouble = paras.mDoubleValues().size();
  size_t check_nMComplex = 0;
  size_t nMComplex = paras.mComplexValues().size();

  // Check size of parameter list
  if (nInt != check_nInt)
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of IntParameters does not match: " +
                       std::to_string(nInt) + " given but " +
                       std::to_string(check_nInt) + " expected."));
  if (nDouble != check_nDouble)
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of FitParameters does not match: " +
                       std::to_string(nDouble) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (nComplex != check_nComplex)
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(nComplex) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (nMInteger != check_nMInteger)
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of MultiInt does not match: " +
                       std::to_string(nMInteger) + " given but " +
                       std::to_string(check_nMInteger) + " expected."));
  if (nMDouble != check_nMDouble)
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of MultiDoubles does not match: " +
                       std::to_string(nMDouble) + " given but " +
                       std::to_string(check_nMDouble) + " expected."));
  if (nMComplex != check_nMComplex)
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of MultiComplexes does not match: " +
                       std::to_string(nMComplex) + " given but " +
                       std::to_string(check_nMComplex) + " expected."));
#endif

  size_t n = paras.mDoubleValue(0)->values().size();
  if (!out)
    out = ComPWA::FunctionTree::MDouble("", n);
  auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
  auto &results = par->values(); // reference
  if (results.size() != n) {
    results.resize(n);
  }

  // Get parameters from ParameterList:
  // We use the same order of the parameters as was used during tree
  // construction.
  unsigned int L = paras.intValue(0)->value();
  double MesonRadius = paras.doubleParameter(0)->value();
  auto s = paras.mDoubleValue(0)->values();
  auto sa = paras.mDoubleValue(1)->values();
  auto sb = paras.mDoubleValue(2)->values();

  if (sa.size() == 1 && sb.size() == 1) {
    double ma = std::sqrt(sa.at(0));
    double mb = std::sqrt(sb.at(0));
    std::transform(paras.mDoubleValue(0)->values().begin(),
                   paras.mDoubleValue(0)->values().end(), results.begin(),
                   [&](double s) {
                     return FormFactorFunctor->operator()(qSquared(s, ma, mb),
                                                          L, MesonRadius);
                   });
  } else if (sa.size() == 1) {
    double ma = std::sqrt(sa.at(0));
    for (size_t i = 0; i < s.size(); ++i) {
      results[i] = FormFactorFunctor->operator()(
          qSquared(s[i], ma, std::sqrt(sb[i])), L, MesonRadius);
    }
  } else if (sb.size() == 1) {
    double mb = std::sqrt(sb.at(0));
    for (size_t i = 0; i < s.size(); ++i) {
      results[i] = FormFactorFunctor->operator()(
          qSquared(s[i], std::sqrt(sa[i]), mb), L, MesonRadius);
    }
  } else {
    for (size_t i = 0; i < s.size(); ++i) {
      results[i] = FormFactorFunctor->operator()(
          qSquared(s[i], std::sqrt(sa[i]), std::sqrt(sb[i])), L, MesonRadius);
    }
  }
}

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA
