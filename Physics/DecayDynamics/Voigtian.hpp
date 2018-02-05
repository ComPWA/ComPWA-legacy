// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef VOIGT_FUNCTION_HPP 
#define VOIGT_FUNCTION_HPP

#include <vector>
#include <memory>
#include <boost/property_tree/ptree.hpp>

#include "Core/Spin.hpp"
#include "Core/Functions.hpp"
#include "Core/Exceptions.hpp"
#include "Physics/DecayDynamics/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"
#include "Physics/DecayDynamics/Utils/Faddeeva.hh"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

class HelicityDecay;

/// \class Voigtian
/// the convolution of a BreitWigner with a Gaussian
/// ref: https://en.wikipedia.org/wiki/Voigtian_profile
///      Voig(x; sigma, gamma) = \int Gaus(x';\sigma)BW(x - x';gamma) dx'
///                            = Re[w(z)]/(\sigma\sqrt{2\pi})
///                              and z = (x + i\gamma)/(\sigma\sqrt{s})
///    : http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
///      this page is a package for computation of w(z)
class Voigtian
    : public ComPWA::Physics::DecayDynamics::AbstractDynamicalFunction {

public:
  Voigtian(std::string name, std::pair<std::string,std::string> daughters,
               std::shared_ptr<ComPWA::PartList> partL);

  /// Check of parameters have changed and normalization has to be recalculatecd
  virtual bool isModified() const;
   
  //================ EVALUATION =================
  std::complex<double> evaluate(const ComPWA::DataPoint &point, int pos) const;

  /// Dynamical Breit-Wigner function.
  /// \param mSq Invariant mass squared
  /// \param mR Mass of the resonant state
  /// \param wR Width of the resonant state
  /// \param sigma Width of the gaussian, i.e., the resolution of the mass spectrum at mR
  /// \return Amplitude value
  static std::complex<double>
  dynamicalFunction(double mSq, double mR, double ma, double mb, double wR,
      unsigned int J, double sigma);

  //============ SET/GET =================

  void SetMassParameter(std::shared_ptr<ComPWA::FitParameter> m) {
    Mass = m;
  }

  std::shared_ptr<ComPWA::FitParameter> GetMassParameter() {
    return Mass;
  }

  void SetMass(double m) { Mass->setValue(m); }

  double GetMass() const { return Mass->value(); }

  void SetWidthParameter(std::shared_ptr<ComPWA::FitParameter> w) {
    Width = w;
  }

  std::shared_ptr<ComPWA::FitParameter> GetWidthParameter() {
    return Width;
  }

  void SetWidth(double w) {Width->setValue(w); }

  double GetWidth() const { return Width->value(); }

  void SetSigma(double sigma) { Sigma = sigma; }

  double GetSigma() const { return Sigma; }

  virtual void parameters(ComPWA::ParameterList &list);

  virtual void parametersFast(std::vector<double> &list) const {
    AbstractDynamicalFunction::parametersFast(list);
    list.push_back(GetMass());
    list.push_back(GetWidth());
  }

  /// Update parameters to the values given in \p par
  virtual void updateParameters(const ComPWA::ParameterList &par);

  //=========== FUNCTIONTREE =================

  virtual bool hasTree() const { return true; }

  virtual std::shared_ptr<ComPWA::FunctionTree>
  tree(const ParameterList &sample, int pos, std::string suffix = "");

protected:
  /// Mass of resonante state
  std::shared_ptr<ComPWA::FitParameter> Mass;
  /// Decay width of resonante state
  std::shared_ptr<ComPWA::FitParameter> Width;
  /// resolution: the width of gaussian function which is used to represent the resolution of mass spectrum
  double Sigma;

private:
  /// Temporary values (used to trigger recalculation of normalization)
  double CurrentMass;
  double CurrentWidth;
};

class VoigtianStrategy : public ComPWA::Strategy {
public:
  VoigtianStrategy(std::string sname = "")
    : ComPWA::Strategy(ParType::MCOMPLEX), name(sname) {}
  
  virtual const std::string to_str() const {
    return ("Voigtian Function of " + name);
  }

  virtual void execute(ComPWA::ParameterList &paras,
                       std::shared_ptr<ComPWA::Parameter> &out); 

protected:
  std::string name;
};

} // namespace DecayDynamics
} // namespace Physics
} // namespace ComPWA

#endif
