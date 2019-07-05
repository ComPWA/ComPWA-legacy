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

#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <vector>

#include "AbstractDynamicalFunction.hpp"
#include "Core/Exceptions.hpp"
#include "Core/FunctionTree/Functions.hpp"
#include "Core/Spin.hpp"
#include "FormFactor.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"
#include "Utils/Faddeeva.hh"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

///
/// \class Voigtian
/// Voigtian class calculate the convolution of a non-relativisitc Breit-Wigner
/// with a Gaussian.
/// ref: https://en.wikipedia.org/wiki/Voigt_profile
///      Voig(x; sigma, gamma) = \int Gaus(x';\sigma)BW(x - x';gamma) dx'
///                            = Re[w(z)]/(\sigma\sqrt{2\pi})
///                              and z = (x + i\gamma)/(\sigma\sqrt{s})
/// In the calculation of voigt function, a Faddeeva Package is used to
/// calculate w(z). ref:
/// http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
///      this page is a package for computation of w(z)
///
class Voigtian : public ComPWA::Physics::Dynamics::AbstractDynamicalFunction {

public:
  //============ CONSTRUCTION ==================
  Voigtian(std::string name, std::pair<std::string, std::string> daughters,
           std::shared_ptr<ComPWA::PartList> partL);

  //================ EVALUATION =================

  std::complex<double> evaluate(const ComPWA::DataPoint &point,
                                unsigned int pos) const;

  /// Dynamical voigt function.
  /// \param mSq Invariant mass squared
  /// \param mR Mass of the resonant state
  /// \param wR Width of the resonant state
  /// \param sigma Width of the gaussian, i.e., the resolution of the mass
  /// spectrum at mR \return Amplitude value
  static std::complex<double> dynamicalFunction(double mSq, double mR,
                                                double wR, double sigma);

  //============ SET/GET =================

  void
  SetWidthParameter(std::shared_ptr<ComPWA::FunctionTree::FitParameter> w) {
    Width = w;
  }

  std::shared_ptr<ComPWA::FunctionTree::FitParameter> GetWidthParameter() {
    return Width;
  }

  void SetWidth(double w) { Width->setValue(w); }

  double GetWidth() const { return Width->value(); }

  void SetMesonRadiusParameter(
      std::shared_ptr<ComPWA::FunctionTree::FitParameter> r) {
    MesonRadius = r;
  }

  std::shared_ptr<ComPWA::FunctionTree::FitParameter>
  GetMesonRadiusParameter() {
    return MesonRadius;
  }

  /// \see GetMesonRadius() const { return MesonRadius->value(); }
  void SetMesonRadius(double w) { MesonRadius->setValue(w); }

  /// Get meson radius.
  /// The meson radius is a measure of the size of the resonant state. It is
  /// used to calculate the angular momentum barrier factors.
  double GetMesonRadius() const { return MesonRadius->value(); }

  /// \see GetFormFactorType()
  void SetFormFactorType(FormFactorType t) { FFType = t; }

  /// Get form factor type.
  /// The type of formfactor that is used to calculate the angular momentum
  /// barrier factors.
  FormFactorType GetFormFactorType() { return FFType; }

  void SetSigma(double sigma) { Sigma = sigma; }

  double GetSigma() const { return Sigma; }

  void updateParametersFrom(const ComPWA::FunctionTree::ParameterList &list);
  void addUniqueParametersTo(ComPWA::FunctionTree::ParameterList &list);
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createFunctionTree(const ComPWA::FunctionTree::ParameterList &DataSample,
                     unsigned int pos, const std::string &suffix) const;

protected:
  /// Orbital Angular Momentum between two daughters in Resonance decay
  ComPWA::Spin L;

  /// Masses of daughter particles
  std::pair<double, double> DaughterMasses;

  /// Names of daughter particles
  std::pair<std::string, std::string> DaughterNames;

  /// Resonance mass
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> Mass;

  /// Decay width of resonante state
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> Width;

  /// Meson radius of resonant state
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> MesonRadius;

  /// Form factor type
  FormFactorType FFType;
  /// resolution: the width of gaussian function which is used to represent the
  /// resolution of mass spectrum
  double Sigma;
};

class VoigtianStrategy : public ComPWA::FunctionTree::Strategy {
public:
  VoigtianStrategy(std::string sname = "")
      : ComPWA::FunctionTree::Strategy(ComPWA::FunctionTree::ParType::MCOMPLEX),
        name(sname) {}

  virtual const std::string to_str() const {
    return ("Voigtian Function of " + name);
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
