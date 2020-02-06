// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// This file contains the declaration of the Voigtian class,
/// which is used the implementation of voigt function, the convolution
/// of a non-relativistic Breit-Wigner and a gaussian.
/// See https://en.wikipedia.org/wiki/voigt_profile for Vogit function.
///

#ifndef VOIGT_FUNCTION_HPP
#define VOIGT_FUNCTION_HPP

#include "Core/FunctionTree/Functions.hpp"
#include "FormFactor.hpp"
#include "RelativisticBreitWigner.hpp"
#include "Utils/Faddeeva.hh"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

///
/// \namespace Voigtian
/// Voigtian is the convolution of a non-relativisitc Breit-Wigner
/// with a Gaussian, see
/// [Wikipedia](https://en.wikipedia.org/wiki/Voigt_profile) \f[
///      \mathrm{Voig}(x; \sigma, \gamma) = \int Gaus(x';\sigma)BW(x - x';gamma)
///      dx'
///                            = Re[w(z)]/(\sigma\sqrt{2\pi})
///                              and z = (x + i\gamma)/(\sigma\sqrt{s})
/// \f]
/// In the calculation of voigt function, a [Faddeeva
/// Package](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package) this page
/// is a package for computation of w(z)) is used to calculate w(z).
///
namespace Voigtian {

struct InputInfo : RelativisticBreitWigner::InputInfo {
  /// resolution: the width of gaussian function which is used to represent the
  /// resolution of mass spectrum
  double Sigma;
};

/// Dynamical voigt function.
/// \param mSq Invariant mass squared
/// \param mR Mass of the resonant state
/// \param wR Width of the resonant state
/// \param sigma Width of the gaussian, i.e., the resolution of the mass
/// spectrum at mR \return Amplitude value
inline std::complex<double> dynamicalFunction(double mSq, double mR, double wR,
                                              double sigma) {

  double sqrtS = sqrt(mSq);

  // the non-relativistic BreitWigner which is convoluted in Voigtian
  // has the exactly following expression:
  // BW(x, m, width) = 1/pi * width/2 * 1/((x - m)^2 + (width/2)^2)
  // i.e., the Lorentz formula with Gamma = width/2 and x' = x - m
  /// https://root.cern.ch/doc/master/RooVoigtianian_8cxx_source.html
  double argu = sqrtS - mR;
  double c = 1.0 / (sqrt(2.0) * sigma);
  double a = c * 0.5 * wR;
  double u = c * argu;
  std::complex<double> z(u, a);
  std::complex<double> v = Faddeeva::w(z, 1e-13);
  double val = c * 1.0 / sqrt(M_PI) * v.real();
  double sqrtVal = sqrt(val);

  /// keep the phi angle of the complex BW
  std::complex<double> invBW(argu, 0.5 * wR);
  std::complex<double> BW = 1.0 / invBW;
  double phi = std::arg(BW);
  std::complex<double> result(sqrtVal * cos(phi), sqrtVal * sin(phi));

  // transform width to coupling
  // Calculate coupling constant to final state
  // MesonRadius = 0.0, noFormFactor
  // std::complex<double> g_final = widthToCoupling(mSq, mR, wR, ma, mb, L, 0.0,
  // formFactorType::noFormFactor);
  // the BW to convolved in voigt is 1/PI * Gamma/2 * 1/((x-m)^2 + (Gamma/2)^2)
  // while I think the one common used in physics is Gamma/2 * 1/((x-m)^2 +
  // (Gamma/2)^2)
  // So we time the PI at last
  std::complex<double> g_final = sqrt(M_PI);
  double g_production = 1;
  result *= g_production;
  result *= g_final;

  assert((!std::isnan(result.real()) || !std::isinf(result.real())) &&
         "Voigtian::dynamicalFunction() | Result is NaN or Inf!");
  assert((!std::isnan(result.imag()) || !std::isinf(result.imag())) &&
         "Voigtian::dynamicalFunction() | Result is NaN or Inf!");

  return result;
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode> createFunctionTree(
    InputInfo Params,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquared);

} // namespace Voigtian

class VoigtianStrategy : public ComPWA::FunctionTree::Strategy {
public:
  VoigtianStrategy(std::string name = "")
      : ComPWA::FunctionTree::Strategy(ComPWA::FunctionTree::ParType::MCOMPLEX,
                                       "Voigtian" + name) {}

  virtual void execute(ComPWA::FunctionTree::ParameterList &paras,
                       std::shared_ptr<ComPWA::FunctionTree::Parameter> &out);
};

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
