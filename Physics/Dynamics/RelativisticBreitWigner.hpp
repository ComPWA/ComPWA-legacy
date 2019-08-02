// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_DYNAMICS_RELATIVISTICBREITWIGNER_HPP_
#define PHYSICS_DYNAMICS_RELATIVISTICBREITWIGNER_HPP_

#include "Core/FunctionTree/FunctionTree.hpp"
#include "Coupling.hpp"
#include "FormFactor.hpp"

#include <vector>

namespace ComPWA {
namespace Physics {
namespace Dynamics {

namespace RelativisticBreitWigner {

struct InputInfo {
  /// Orbital Angular Momentum between two daughters in Resonance decay
  unsigned int L;
  /// Masses of daughter particles
  std::pair<std::shared_ptr<ComPWA::FunctionTree::FitParameter>,
            std::shared_ptr<ComPWA::FunctionTree::FitParameter>>
      DaughterMasses;
  /// Resonance mass
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> Mass;
  /// Decay width of resonant state
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> Width;
  /// Meson radius of resonant state
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> MesonRadius;
  /// Form factor type
  FormFactorType FFType;
};

///
/// Relativistic Breit-Wigner model with barrier factors.
/// The dynamical function implemented here is taken from PDG2014 (Eq.47.22)
/// for the one channel case. The dynamic reaction
/// \f[
/// \mathcal{A}_R(s) = \frac{g_p*g}{s - M_R^2 + i \sqrt{s} \Gamma_R B^2}
/// \f]
/// \f$ g_p, g\f$ are the coupling constants for production and decay and
/// the barrier term \f$ B^2\f$ is parameterized according to Eq.47.23:
/// \f[
///     B^2 = \left( \frac{q(s)}{q(s_R)} \right)^{2L+1} \times
//                     \left( \frac{F(s)}{F(s_R)} \right)^{2}
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
/// \param ffType Form factor type
/// \return Amplitude value
inline std::complex<double>
dynamicalFunction(double mSq, double mR, double ma, double mb, double width,
                  unsigned int L, double mesonRadius, FormFactorType ffType) {

  std::complex<double> i(0, 1);
  double sqrtS = std::sqrt(mSq);

  // Phase space factors at sqrt(s) and at the resonance position
  auto phspFactorSqrtS = phspFactor(sqrtS, ma, mb);
  auto phspFactormR = phspFactor(mR, ma, mb);

  // Check if we have an event which is exactly at the phase space boundary
  if (phspFactorSqrtS == std::complex<double>(0, 0))
    return std::complex<double>(0, 0);

  std::complex<double> qRatio =
      std::pow((phspFactorSqrtS / phspFactormR) * mR / sqrtS, (2 * L + 1));
  double ffR = FormFactor(mR, ma, mb, L, mesonRadius, ffType);
  double ff = FormFactor(sqrtS, ma, mb, L, mesonRadius, ffType);
  std::complex<double> barrierTermSq = qRatio * (ff * ff) / (ffR * ffR);

  // Calculate normalized vertex function gammaA(s_R) at the resonance position
  // (see PDG2014, Chapter 47.2)
  std::complex<double> gammaA(1, 0); // spin==0
  if (L > 0) {
    std::complex<double> qR = std::pow(qValue(mR, ma, mb), L);
    gammaA = ffR * qR;
  }

  // Coupling to the final state (ma, mb)
  std::complex<double> g_final =
      widthToCoupling(mR, width, gammaA, phspFactorSqrtS);

  std::complex<double> denom(mR * mR - mSq, 0);
  denom += (-1.0) * i * sqrtS * (width * barrierTermSq);

  std::complex<double> result = g_final / denom;

  assert(
      (!std::isnan(result.real()) || !std::isinf(result.real())) &&
      "RelativisticBreitWigner::dynamicalFunction() | Result is NaN or Inf!");
  assert(
      (!std::isnan(result.imag()) || !std::isinf(result.imag())) &&
      "RelativisticBreitWigner::dynamicalFunction() | Result is NaN or Inf!");

  return result;
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
createFunctionTree(InputInfo Params,
                   const ComPWA::FunctionTree::ParameterList &DataSample,
                   unsigned int pos, std::string suffix);

} // namespace RelativisticBreitWigner

class BreitWignerStrategy : public ComPWA::FunctionTree::Strategy {
public:
  BreitWignerStrategy(std::string namee = "")
      : ComPWA::FunctionTree::Strategy(ComPWA::FunctionTree::ParType::MCOMPLEX),
        name(namee) {}

  virtual const std::string to_str() const {
    return ("relativistic BreitWigner of " + name);
  }

  virtual void execute(ComPWA::FunctionTree::ParameterList &paras,
                       std::shared_ptr<ComPWA::FunctionTree::Parameter> &out);

protected:
  std::string name;
};

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
