// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Flatte.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

using ComPWA::FunctionTree::FitParameter;
using ComPWA::FunctionTree::Parameter;
using ComPWA::FunctionTree::ParameterList;
using ComPWA::FunctionTree::TreeNode;
using ComPWA::FunctionTree::Value;

std::shared_ptr<TreeNode> Flatte::createFunctionTree(
    InputInfo Params,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquared) {
  if (Params.HiddenCouplings.size() == 1)
    Params.HiddenCouplings.push_back(Coupling(0.0, 0.0, 0.0));

  if (Params.HiddenCouplings.size() != 2)
    throw std::runtime_error(
        "Flatte::createFunctionTree() | Vector with "
        "couplings has a wrong size. We expect either 2 or 3 couplings.");

  if (!Params.G)
    throw std::runtime_error(
        "Flatte::createFunctionTree() | Coupling to signal channel not set");

  size_t sampleSize = InvMassSquared->values().size();

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<TreeNode>(
      MComplex("", sampleSize),
      std::make_shared<FlatteStrategy>(Params.FormFactorFunctor));

  tr->addNodes({createLeaf(Params.Mass), createLeaf(Params.G)});
  for (unsigned int i = 0; i < Params.HiddenCouplings.size(); ++i) {
    tr->addNodes({createLeaf(Params.HiddenCouplings.at(i).MassA),
                  createLeaf(Params.HiddenCouplings.at(i).MassB),
                  createLeaf(Params.HiddenCouplings.at(i).G)});
  }
  tr->addNodes({createLeaf((double)Params.L, "L"),
                createLeaf(Params.MesonRadius), createLeaf(InvMassSquared),
                createLeaf(Params.DaughterInvariantMasses.first),
                createLeaf(Params.DaughterInvariantMasses.second)});

  return tr;
}

void FlatteStrategy::execute(ParameterList &paras,
                             std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("FlatteStrategy::execute() | Parameter type mismatch!");

#ifndef NDEBUG
  // Check parameter type
  if (checkType != out->type())
    throw(
        WrongParType("FlatteStrategy::execute() | "
                     "Output parameter is of type " +
                     std::string(ComPWA::FunctionTree::ParNames[out->type()]) +
                     " and conflicts with expected type " +
                     std::string(ComPWA::FunctionTree::ParNames[checkType])));

  // How many parameters do we expect?
  size_t check_nInt = 0;
  size_t nInt = paras.intValues().size();
  size_t check_nDouble = 10;
  size_t nDouble = paras.doubleValues().size();
  nDouble += paras.doubleParameters().size();
  size_t check_nComplex = 0;
  size_t nComplex = paras.complexValues().size();
  size_t check_nMInteger = 0;
  size_t nMInteger = paras.mIntValues().size();
  size_t check_nMDouble = 3;
  size_t nMDouble = paras.mDoubleValues().size();
  size_t check_nMComplex = 0;
  size_t nMComplex = paras.mComplexValues().size();

  // Check size of parameter list
  if (nInt != check_nInt)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of IntParameters does not match: " +
                       std::to_string(nInt) + " given but " +
                       std::to_string(check_nInt) + " expected."));
  if (nDouble != check_nDouble)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of FitParameters does not match: " +
                       std::to_string(nDouble) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (nComplex != check_nComplex)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(nComplex) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (nMInteger != check_nMInteger)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of MultiInt does not match: " +
                       std::to_string(nMInteger) + " given but " +
                       std::to_string(check_nMInteger) + " expected."));
  if (nMDouble != check_nMDouble)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of MultiDoubles does not match: " +
                       std::to_string(nMDouble) + " given but " +
                       std::to_string(check_nMDouble) + " expected."));
  if (nMComplex != check_nMComplex)
    throw(BadParameter("FlatteStrategy::execute() | "
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
  // calc function for each point
  for (size_t ele = 0; ele < n; ele++) {
    try {
      // Generally we need to add a factor q^{2J+1} to each channel term.
      // But since Flatte resonances are usually J=0 we neglect it here.
      results.at(ele) = Flatte::dynamicalFunction(
          paras.mDoubleValue(0)->values().at(ele),
          paras.doubleParameter(0)->value(),  // mass
          paras.doubleParameter(1)->value(),  // g1_massA
          paras.doubleParameter(2)->value(),  // g1_massB
          paras.doubleParameter(3)->value(),  // g1
          paras.doubleParameter(4)->value(),  // g2_massA
          paras.doubleParameter(5)->value(),  // g2_massB
          paras.doubleParameter(6)->value(),  // g2
          paras.doubleParameter(7)->value(),  // g3_massA
          paras.doubleParameter(8)->value(),  // g3_massB
          paras.doubleParameter(9)->value(),  // g3
          paras.doubleValue(0)->value(),      // OrbitalAngularMomentum
          paras.doubleParameter(10)->value(), // mesonRadius
          FormFactorFunctor);
    } catch (std::exception &ex) {
      LOG(ERROR) << "FlatteStrategy::execute() | " << ex.what();
      throw(std::runtime_error("FlatteStrategy::execute() | "
                               "Evaluation of dynamic function failed!"));
    }
  }
}

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA
