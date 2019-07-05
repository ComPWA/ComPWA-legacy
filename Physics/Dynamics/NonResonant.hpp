// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_DYNAMICS_NONRESONANT_HPP_
#define COMPWA_PHYSICS_DYNAMICS_NONRESONANT_HPP_

#include "AbstractDynamicalFunction.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

class NonResonant : public AbstractDynamicalFunction {

public:
  NonResonant(std::string name = "") : AbstractDynamicalFunction(name){};
  virtual ~NonResonant(){};

  /*boost::property_tree::ptree save() const {
    return boost::property_tree::ptree();
  }*/

  std::complex<double> evaluate(const ComPWA::DataPoint &point,
                                unsigned int pos) const;

  void updateParametersFrom(const ComPWA::FunctionTree::ParameterList &list);
  void addUniqueParametersTo(ComPWA::FunctionTree::ParameterList &list);
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createFunctionTree(const ComPWA::FunctionTree::ParameterList &DataSample,
                     unsigned int pos, const std::string &suffix) const;
};

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
