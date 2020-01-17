// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef Coupling_h
#define Coupling_h

#include "Core/FunctionTree/FitParameter.hpp"
#include "Core/Properties.hpp"
#include "FormFactor.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

/// Convert width to complex coupling. This is the implementation of
/// PDG2014, Chapter 47.2, Eq. 47.21 (inverted).
inline std::complex<double> couplingToWidth(double mR, double g,
                                            std::complex<double> gamma,
                                            std::complex<double> phspFactor) {
  // calculate phsp factor
  std::complex<double> res = std::norm(gamma) * g * g * phspFactor / mR;

#ifndef NDEBUG
  // check for NaN
  if (std::isnan(res.real()) || std::isnan(res.imag()))
    throw std::runtime_error("couplingToWidth() | Result is NaN!");
  // check for inf
  if (std::isinf(res.real()) || std::isinf(res.imag()))
    throw std::runtime_error("couplingToWidth() | Result is inf!");
#endif

  return res;
}

/// Convert width to complex coupling. The form factor \p formFactorR, the
/// normalized vertex function \p gamma (both evaluated at the resonance pole)
/// and the \p phspFactor (evaluated at sqrt(s)) can be passed in
/// order to save computation time. This is the implementation of
/// PDG2014, Chapter 47.2, Eq. 47.21. See also widthToCoupling(double mSq,
/// double mR, double width, double ma, double mb, double spin, double
/// mesonRadius, formFactorType type)
inline std::complex<double> widthToCoupling(double mR, double width,
                                            std::complex<double> gamma,
                                            std::complex<double> phspFactor) {

  auto denom = gamma * std::sqrt(phspFactor);
  auto res = std::complex<double>(sqrt(mR * width), 0) / denom;

#ifndef NDEBUG
  // check for NaN
  if (std::isnan(res.real()) || std::isnan(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::widthToCoupling() | "
                             "Result is NaN!");
  // check for inf
  if (std::isinf(res.real()) || std::isinf(res.imag()))
    throw std::runtime_error("AbstractDynamicalFunction::widthToCoupling() | "
                             "Result is inf!");
#endif
  return res;
}

class Coupling {
public:
  Coupling(){};
  Coupling(double c, double massA, double massB)
      : G(new ComPWA::FunctionTree::FitParameter("", c)),
        MassA(new ComPWA::FunctionTree::FitParameter("", massA)),
        MassB(new ComPWA::FunctionTree::FitParameter("", massB)){};

  Coupling(const ComPWA::ParticleList &partL,
           const boost::property_tree::ptree tr) {
    G = std::make_shared<ComPWA::FunctionTree::FitParameter>(tr.get_child(""));
    std::string nameA = tr.get<std::string>("ParticleA");
    std::string nameB = tr.get<std::string>("ParticleB");
    auto mA = findParticle(partL, nameA).getMass();
    auto mB = findParticle(partL, nameB).getMass();
    MassA = std::make_shared<ComPWA::FunctionTree::FitParameter>(mA);
    MassB = std::make_shared<ComPWA::FunctionTree::FitParameter>(mB);
  };

  std::shared_ptr<ComPWA::FunctionTree::FitParameter> G;
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> MassA;
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> MassB;
};

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
