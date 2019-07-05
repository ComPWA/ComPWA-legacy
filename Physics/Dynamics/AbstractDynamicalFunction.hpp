// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_DYNAMICS_ABSTRACTDYNAMICALFUNCTION_HPP_
#define COMPWA_PHYSICS_DYNAMICS_ABSTRACTDYNAMICALFUNCTION_HPP_

#include <complex>

#include <boost/property_tree/ptree.hpp>

#include "Core/Event.hpp"
#include "Core/FunctionTree/FunctionTree.hpp"
#include "Core/FunctionTree/ParameterList.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

class AbstractDynamicalFunction : public ComPWA::FunctionTree::Optimizable {

public:
  AbstractDynamicalFunction(std::string name = "") : Name(name){};

  virtual ~AbstractDynamicalFunction(){};

  virtual std::complex<double> evaluate(const ComPWA::DataPoint &point,
                                        unsigned int pos) const = 0;

  virtual std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createFunctionTree(const ComPWA::FunctionTree::ParameterList &DataSample,
                     unsigned int pos, const std::string &suffix) const = 0;

protected:
  std::string Name;
};

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
