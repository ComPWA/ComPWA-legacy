// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_DYNAMICS_RELATIVISTICBREITWIGNER_HPP_
#define PHYSICS_DYNAMICS_RELATIVISTICBREITWIGNER_HPP_

#include "Core/FunctionTree/TreeNode.hpp"
#include "Coupling.hpp"
#include "FormFactor.hpp"

#include <vector>

namespace ComPWA {
namespace Physics {
namespace Dynamics {

struct InputInfo {
  std::string Type;
  /// Orbital Angular Momentum between two daughters in Resonance decay
  unsigned int L;
  /// Invariant Masses of daughter particles
  std::pair<std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>,
            std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>>
      DaughterInvariantMasses;
  /// Resonance mass
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> Mass;

  /// Meson radius of resonant state
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> MesonRadius;
  /// Form factor function object
  std::shared_ptr<FormFactor> FormFactorFunctor;
};

namespace RelativisticBreitWigner {

using BreitWignerFunction = std::function<std::complex<double>(
    double, double, double, double, double, unsigned int, double,
    std::shared_ptr<FormFactor>)>;

struct InputInfo : ComPWA::Physics::Dynamics::InputInfo {
  /// Decay width of resonant state
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> Width;
};

///
/// Relativistic Breit-Wigner model with barrier factors.
/// The dynamical function implemented here is taken from PDG2018 (Eq.48.22)
/// for the one channel case. The dynamic reaction
/// \f[
/// \mathcal{A}_R(s) = \frac{g_p*g}{s - M_R^2 + i \sqrt{s} \Gamma_R B^2}
/// \f]
/// \f$ g_p, g\f$ are the coupling constants for production and decay and
/// the barrier term \f$ B^2\f$ is parameterized according to Eq.48.23:
/// \f[
///     B^2 = \left( \frac{q(\sqrt{s})}{q(M_R)} \right)^{2L+1} \times
///           \left( \frac{M_R}{\sqrt{s}} \right) \times
///           \left( \frac{F(\sqrt{s})}{F(\sqrt{s_R})} \right)^{2}
/// \f]
/// This corresponds to the Blatt Weisskopf form factors B_L like
/// \f[
///     B^2 = \left( \frac{q(\sqrt{s})}{q(M_R)} \right) \times
///           \left( \frac{M_R}{\sqrt{s}} \right) \times
///           \left( \frac{B_L(\sqrt{s})}{B_L(\sqrt{s_R})} \right)^{2}
/// \f]
///
/// \param mSq Invariant mass squared
/// \param mR Mass of the resonant state
/// \param ma Mass of daughter particle
/// \param mb Mass of daughter particle
/// \param width Decay width
/// \param L Orbital angular momentum between two daughters a and b
/// \param mesonRadius Meson Radius
/// \param FormFactorFunctor Form factor functor
/// \return Amplitude value
inline std::complex<double>
relativisticBreitWigner(double mSq, double mR, double ma, double mb,
                        double width, unsigned int L, double mesonRadius,
                        std::shared_ptr<FormFactor> FormFactorFunctor) {

  double sqrtS = std::sqrt(mSq);

  // Phase space factors at sqrt(s) and at the resonance position
  auto phspFactorSqrtS = phspFactor(sqrtS, ma, mb);
  auto phspFactormR = phspFactor(mR, ma, mb);

  // Check if we have an event which is exactly at the phase space boundary
  if (phspFactorSqrtS == std::complex<double>(0, 0))
    return std::complex<double>(0, 0);

  // The each FormFactor includes a term q^L. Therefore only q instead of
  // q^(2L+1) has to be calculated. Also because phspFactor ~ q/sqrt(s) (see
  // PDG2018, equation 48.2) the factor \frac{M_R}{\sqrt{s}} can also be omitted
  double qRatio = (phspFactorSqrtS / phspFactormR);

  double ffR =
      FormFactorFunctor->operator()(qSquared(mR * mR, ma, mb), L, mesonRadius);
  double ff =
      FormFactorFunctor->operator()(qSquared(mSq, ma, mb), L, mesonRadius);
  double barrierTermSq = qRatio * (ff * ff) / (ffR * ffR);

  std::complex<double> denom(mR * mR - mSq,
                             -1. * sqrtS * width * barrierTermSq);

  std::complex<double> result = mR * width / denom;

  assert((!std::isnan(result.real()) || !std::isinf(result.real())) &&
         "RelativisticBreitWigner::standardRelBW() | Result is NaN or Inf!");
  assert((!std::isnan(result.imag()) || !std::isinf(result.imag())) &&
         "RelativisticBreitWigner::standardRelBW() | Result is NaN or Inf!");

  return result;
}

///
/// Relativistic Breit-Wigner model with barrier factors.
/// The dynamical function implemented here is taken from PDG2018 (Eq.48.22)
/// for the one channel case. The dynamic reaction
/// \f[
/// \mathcal{A}_R(s) = \frac{g_p*g}{s - M_R^2 + i \sqrt{s} \Gamma_R B^2}
/// \f]
/// \f$ g_p, g\f$ are the coupling constants for production and decay and
/// the barrier term \f$ B^2\f$ is parameterized according to Eq.48.23:
/// \f[
///     B^2 = \left( \frac{q(\sqrt{s})}{q(M_R)} \right)^{2L+1} \times
///           \left( \frac{M_R}{\sqrt{s}} \right) \times
///           \left( \frac{F(\sqrt{s})}{F(\sqrt{s_R})} \right)^{2}
/// \f]
/// This corresponds to the Blatt Weisskopf form factors B_L like
/// \f[
///     B^2 = \left( \frac{q(\sqrt{s})}{q(M_R)} \right) \times
///           \left( \frac{M_R}{\sqrt{s}} \right) \times
///           \left( \frac{B_L(\sqrt{s})}{B_L(\sqrt{s_R})} \right)^{2}
/// \f]
///
/// \param mSq Invariant mass squared
/// \param mR Mass of the resonant state
/// \param ma Mass of daughter particle
/// \param mb Mass of daughter particle
/// \param width Decay width
/// \param L Orbital angular momentum between two daughters a and b
/// \param mesonRadius Meson Radius
/// \param FormFactorFunctor Form factor functor
/// \return Amplitude value
inline std::complex<double> relativisticBreitWignerAnalyticCont(
    double mSq, double mR, double ma, double mb, double width, unsigned int L,
    double mesonRadius, std::shared_ptr<FormFactor> FormFactorFunctor) {

  std::complex<double> i(0, 1);
  double sqrtS = std::sqrt(mSq);

  // Phase space factors at sqrt(s) and at the resonance position
  auto phspFactorSqrtS = phspFactorAC(sqrtS, ma, mb);
  auto phspFactormR = phspFactorAC(mR, ma, mb);

  // Check if we have an event which is exactly at the phase space boundary
  if (phspFactorSqrtS == std::complex<double>(0, 0))
    return std::complex<double>(0, 0);

  // The each FormFactor includes a term q^L. Therefore only q instead of
  // q^(2L+1) has to be calculated. Also because phspFactor ~ q/sqrt(s) (see
  // PDG2018, equation 48.2) the factor \frac{M_R}{\sqrt{s}} can also be omitted
  std::complex<double> qRatio = (phspFactorSqrtS / phspFactormR);

  double ffR = FormFactorFunctor->operator()(
      std::pow(std::abs(qValueAC(mR, ma, mb)), 2), L, mesonRadius);
  double ff = FormFactorFunctor->operator()(
      std::pow(std::abs(qValueAC(sqrtS, ma, mb)), 2), L, mesonRadius);
  std::complex<double> barrierTermSq = qRatio * (ff * ff) / (ffR * ffR);

  // Calculate normalized vertex function gammaA(s_R) at the resonance position
  // (see PDG2018, Chapter 48.2)
  std::complex<double> gammaA(1, 0); // spin==0
  if (L > 0) {
    std::complex<double> qR = std::pow(qValueAC(mR, ma, mb), L);
    gammaA = ffR * qR;
  }

  // Coupling to the final state (ma, mb)
  std::complex<double> g_final =
      widthToCoupling(mR, width, gammaA, phspFactorSqrtS);

  std::complex<double> denom(mR * mR - mSq, 0);
  denom += (-1.0) * i * sqrtS * (width * barrierTermSq);

  std::complex<double> result = g_final / denom;

  assert((!std::isnan(result.real()) || !std::isinf(result.real())) &&
         "RelativisticBreitWigner::relativisticBreitWignerAnalyticCont() | "
         "Result is NaN or Inf!");
  assert((!std::isnan(result.imag()) || !std::isinf(result.imag())) &&
         "RelativisticBreitWigner::relativisticBreitWignerAnalyticCont() | "
         "Result is NaN or Inf!");

  return result;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode> createFunctionTree(
    InputInfo Params,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquared);

} // namespace RelativisticBreitWigner

class BreitWignerStrategy : public ComPWA::FunctionTree::Strategy {
public:
  BreitWignerStrategy(std::shared_ptr<FormFactor> FF,
                      RelativisticBreitWigner::BreitWignerFunction BWFunction_)
      : ComPWA::FunctionTree::Strategy(ComPWA::FunctionTree::ParType::MCOMPLEX,
                                       "RelativisticBreitWigner"),
        FormFactorFunctor(FF), BWFunction(BWFunction_) {}

  virtual void execute(ComPWA::FunctionTree::ParameterList &paras,
                       std::shared_ptr<ComPWA::FunctionTree::Parameter> &out);

private:
  std::shared_ptr<FormFactor> FormFactorFunctor;
  RelativisticBreitWigner::BreitWignerFunction BWFunction;
};

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
