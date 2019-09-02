// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "FormFactor.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

using ComPWA::FunctionTree::FitParameter;
using ComPWA::FunctionTree::FunctionTree;
using ComPWA::FunctionTree::Parameter;
using ComPWA::FunctionTree::ParameterList;
using ComPWA::FunctionTree::Value;

std::shared_ptr<ComPWA::FunctionTree::FunctionTree> createFunctionTree(
    std::string Name,
    std::shared_ptr<ComPWA::FunctionTree::FitParameter> Daughter1Mass,
    std::shared_ptr<ComPWA::FunctionTree::FitParameter> Daughter2Mass,
    std::shared_ptr<ComPWA::FunctionTree::FitParameter> MesonRadius,
    unsigned int L, FormFactorType FFType, const ParameterList &DataSample,
    unsigned int pos, std::string suffix) {
  size_t sampleSize = DataSample.mDoubleValue(0)->values().size();

  std::string ffNodeName = "ProductionFormFactor(" + Name + ")" + suffix;
  auto ffTree = std::make_shared<FunctionTree>(
      ffNodeName, ComPWA::FunctionTree::MDouble("", sampleSize),
      std::make_shared<FormFactorStrategy>());
  // add L and FFType as double value leaf, since there is no int leaf
  ffTree->createLeaf("OrbitalAngularMomentum", L, ffNodeName);
  ffTree->createLeaf("MesonRadius", MesonRadius, ffNodeName);
  ffTree->createLeaf("FormFactorType", (double)FFType, ffNodeName);
  ffTree->createLeaf("MassA", Daughter1Mass, ffNodeName);
  ffTree->createLeaf("MassB", Daughter2Mass, ffNodeName);
  ffTree->createLeaf("Data_mSq[" + std::to_string(pos) + "]",
                     DataSample.mDoubleValue(pos), ffNodeName);
  ffTree->parameter();

  if (!ffTree->sanityCheck())
    throw std::runtime_error("ProductionFormFactor::createFunctionTree | "
                             "Tree didn't pass sanity check!");

  return ffTree;
};

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
  size_t check_nInt = 0;
  size_t nInt = paras.intValues().size();
  // L, MesonRadius, FFType, Daughter1Mass, Daughter2Mass
  size_t check_nDouble = 5;
  size_t nDouble = paras.doubleValues().size();
  nDouble += paras.doubleParameters().size();
  size_t check_nComplex = 0;
  size_t nComplex = paras.complexValues().size();
  size_t check_nMInteger = 0;
  size_t nMInteger = paras.mIntValues().size();
  // DataSample.mDoubleValue(pos) (mSq)
  size_t check_nMDouble = 1;
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
  unsigned int orbitL = paras.doubleValue(0)->value();
  double MesonRadius = paras.doubleParameter(0)->value();
  FormFactorType ffType = FormFactorType(paras.doubleValue(1)->value());
  double ma = paras.doubleParameter(1)->value();
  double mb = paras.doubleParameter(2)->value();

  // calc function for each point
  for (unsigned int ele = 0; ele < n; ele++) {
    try {
      results.at(ele) = FormFactor(paras.mDoubleValue(0)->values().at(ele), ma,
                                   mb, orbitL, MesonRadius, ffType);
    } catch (std::exception &ex) {
      LOG(ERROR) << "FormFactorStrategy::execute() | " << ex.what();
      throw(std::runtime_error("FormFactorStrategy::execute() | "
                               "Evaluation of dynamic function failed!"));
    }
  }
}

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA
