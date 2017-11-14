

// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Estimator/MinLogLH/SumMinLogLH.hpp"

using namespace ComPWA::Estimator;

SumMinLogLH::SumMinLogLH() : _nCalls(0) {}

double SumMinLogLH::ControlParameter(ParameterList &minPar) {
  double lh = 0;
  if (!_tree) {
    for (auto i : _minLogLh)
      lh = +i->ControlParameter(minPar);
  } else {
    _tree->Recalculate();
    auto logLH =
        std::dynamic_pointer_cast<DoubleParameter>(_tree->Head()->Parameter());
    lh = logLH->GetValue();
  }
  _nCalls++;
  return lh; // return -logLH
}

void SumMinLogLH::UseFunctionTree(bool onoff) {
  if (onoff && _tree)
    return;     // Tree already exists
  if (!onoff) { // disable tree
    _tree = std::shared_ptr<FunctionTree>();
    return;
  }
  _tree = std::make_shared<FunctionTree>(
      "SumLogLh", std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)));
  for (auto tr : _minLogLh) {
    try {
      tr->UseFunctionTree(true);
    } catch (std::exception &ex) {
      LOG(error)
          << "SumMinLogLH::UseFunctionTree() | Construction of one or more sub "
             "trees has failed! Error: "
          << ex.what();
      throw;
    }
    _tree->InsertTree(tr->GetTree(), "SumLogLh");
  }

  _tree->Recalculate();
  if (!_tree->SanityCheck()) {
    throw std::runtime_error(
        "SumMinLogLH::UseFunctionTree() | Tree has structural "
        "problems. Sanity check not passed!");
  }
  return;
}

std::shared_ptr<ComPWA::FunctionTree> SumMinLogLH::GetTree() { return _tree; }
