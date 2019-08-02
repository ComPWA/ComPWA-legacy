// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Flatte.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

using ComPWA::FunctionTree::FitParameter;
using ComPWA::FunctionTree::FunctionTree;
using ComPWA::FunctionTree::Parameter;
using ComPWA::FunctionTree::ParameterList;
using ComPWA::FunctionTree::Value;

std::shared_ptr<FunctionTree> Flatte::createFunctionTree(
    InputInfo Params, const ComPWA::FunctionTree::ParameterList &DataSample,
    unsigned int pos, std::string suffix) {

  auto Couplings = Params.Couplings;

  if (Couplings.size() != 2 && Couplings.size() != 3)
    throw std::runtime_error(
        "AmpFlatteRes::SetCouplings() | Vector with "
        "couplings has a wrong size. We expect either 2 or 3 couplings.");

  if (Couplings.size() == 2)
    Couplings.push_back(Coupling(0.0, 0.0, 0.0));
  // Check if one of the  coupling match the final state (_daughterMasses)
  if (!(Params.DaughterMasses.first && Params.DaughterMasses.second))
    LOG(INFO)
        << "AmpFlatteRes::SetCouplings() | Masses of decay products not set. "
           " Can not determine if correct couplings were set.";

  bool ok = false;
  for (auto i : Couplings) {
    if (i.GetMassA() == Params.DaughterMasses.first->value() &&
        i.GetMassB() == Params.DaughterMasses.second->value())
      ok = true;
    if (i.GetMassB() == Params.DaughterMasses.first->value() &&
        i.GetMassA() == Params.DaughterMasses.second->value())
      ok = true;
  }
  if (!ok)
    throw std::runtime_error("AmpFlatteRes::SetCouplings() | No couplings "
                             "for the current decay particles set!");

  size_t sampleSize = DataSample.mDoubleValue(pos)->values().size();

  std::string NodeName = "Flatte" + suffix;
  auto tr = std::make_shared<FunctionTree>(
      NodeName, ComPWA::FunctionTree::MComplex("", sampleSize),
      std::make_shared<FlatteStrategy>(""));

  tr->createLeaf("Mass", Params.Mass, NodeName);
  for (unsigned int i = 0; i < Params.Couplings.size(); ++i) {
    tr->createLeaf("g_" + std::to_string(i) + "_massA",
                   Params.Couplings.at(i).GetMassA(), NodeName);
    tr->createLeaf("g_" + std::to_string(i) + "_massB",
                   Params.Couplings.at(i).GetMassB(), NodeName);
    tr->createLeaf("g_" + std::to_string(i),
                   Params.Couplings.at(i).GetValueParameter(), NodeName);
  }
  tr->createLeaf("OrbitalAngularMomentum", Params.L, NodeName);
  tr->createLeaf("MesonRadius", Params.MesonRadius, NodeName);
  tr->createLeaf("FormFactorType", Params.FFType, NodeName);
  tr->createLeaf(DataSample.mDoubleValue(pos)->name(),
                 DataSample.mDoubleValue(pos), NodeName);

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
  size_t check_nDouble = 13;
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
          paras.doubleParameter(0)->value(), // mass
          paras.doubleValue(0)->value(),     // g1_massA
          paras.doubleValue(1)->value(),     // g1_massB
          paras.doubleParameter(1)->value(), // g1
          paras.doubleValue(2)->value(),     // g2_massA
          paras.doubleValue(3)->value(),     // g2_massB
          paras.doubleParameter(2)->value(), // g2
          paras.doubleValue(4)->value(),     // g3_massA
          paras.doubleValue(5)->value(),     // g3_massB
          paras.doubleParameter(3)->value(), // g3
          paras.doubleValue(6)->value(),     // OrbitalAngularMomentum
          paras.doubleParameter(4)->value(), // mesonRadius
          FormFactorType(paras.doubleValue(7)->value()) // ffType
      );
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
