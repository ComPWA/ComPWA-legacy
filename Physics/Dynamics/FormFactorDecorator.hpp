// Copyright (c) 2013, 2019 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef FORMFACTOR_DECORATOR_HPP_
#define FORMFACTOR_DECORATOR_HPP_

#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <vector>

#include "AbstractDynamicalFunction.hpp"
#include "Core/Exceptions.hpp"
#include "Core/FunctionTree/Functions.hpp"
#include "Core/Spin.hpp"
#include "FormFactor.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

/// The production formfactor is implemented in Blatt-Weisskopf type
/// \f$B_{L}(q)\f$ Ref. Dalitz Plot Analysis, PDG 2012, where \f[
/// \frac{B^{\prime}(q_{0})}{B^{\prime}_{L}(q)} = q^{L}B^{\prime}_{L}(q, q_{0})
/// \f]
/// \f$B_{L}(q)(i.e., F_L(q))\f$ with \f$L\f$ up to 4 can be found in
/// BNL-QGS-06-101
class FormFactorDecorator : public AbstractDynamicalFunction {

public:
  //============ CONSTRUCTION ==================
  FormFactorDecorator(
      std::string name,
      std::shared_ptr<AbstractDynamicalFunction> undecoratedBreitWigner,
      std::shared_ptr<ComPWA::FunctionTree::FitParameter> mass1,
      std::shared_ptr<ComPWA::FunctionTree::FitParameter> mass2,
      std::shared_ptr<ComPWA::FunctionTree::FitParameter> radius,
      ComPWA::Spin orbitL, FormFactorType ffType);
  virtual ~FormFactorDecorator();

  //================ EVALUATION =================

  std::complex<double> evaluate(const ComPWA::DataPoint &point,
                                unsigned int pos) const;

  /// Blatt-Weisskopf formfactors for production of R -> a b.
  /// \param mSq Invariant mass squared
  /// \param ma Mass of daughter particle
  /// \param mb Mass of daughter particle
  /// \param L Orbital angular momentum between two daughters a and b
  /// \param mesonRadius Meson Radius
  /// \param ffType Form factor type
  static double formFactor(double mSq, double ma, double mb, unsigned int L,
                           double mesonRadius, FormFactorType ffType);

  void updateParametersFrom(const ComPWA::FunctionTree::ParameterList &list);
  void addUniqueParametersTo(ComPWA::FunctionTree::ParameterList &list);
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createFunctionTree(const ComPWA::FunctionTree::ParameterList &DataSample,
                     unsigned int pos, const std::string &suffix) const;

private:
  std::string Name;
  std::shared_ptr<AbstractDynamicalFunction> UndecoratedBreitWigner;
  /// Mass of daughters
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> Daughter1Mass;
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> Daughter2Mass;
  /// Meson radius of resonant state
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> MesonRadius;

  /// Orbital Angular Momentum between two daughters in Resonance decay
  ComPWA::Spin L;
  /// Form factor type
  FormFactorType FFType;
};

class FormFactorStrategy : public ComPWA::FunctionTree::Strategy {
public:
  FormFactorStrategy(std::string namee = "")
      : ComPWA::FunctionTree::Strategy(ComPWA::FunctionTree::ParType::MDOUBLE),
        name(namee) {}

  virtual const std::string to_str() const {
    return ("production FormFactorStratey of " + name);
  }

  virtual void execute(ComPWA::FunctionTree::ParameterList &paras,
                       std::shared_ptr<ComPWA::FunctionTree::Parameter> &out);

private:
  std::string name;
};

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
