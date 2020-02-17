// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_DYNAMICS_FORMFACTOR_HPP_
#define COMPWA_PHYSICS_DYNAMICS_FORMFACTOR_HPP_

#include <complex>
#include <exception>

#include "Core/FunctionTree/Functions.hpp"
#include "Core/FunctionTree/TreeNode.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

/// Calculate Break-up momentum squared.
/// At energy \p S for particles with masses \p sqrtSA and \p sqrtSB. From
/// PDG2014 Eq.46-20a.
/// \param S squared invariant mass of decaying system
/// \param sqrtSA invariant mass of decay product A
/// \param sqrtSB invariant mass of decay product B
inline double qSquared(double S, double sqrtSA, double sqrtSB) {
  double mapb = sqrtSA + sqrtSB;
  double mamb = sqrtSA - sqrtSB;
  double t1 = S - mapb * mapb;
  double t2 = S - mamb * mamb;
  return (t1 * t2 / (4 * S));
}

inline double phspFactor(double sqrtS, double ma, double mb) {
  return std::abs(std::sqrt(qSquared(sqrtS * sqrtS, ma, mb))) /
         (8.0 * M_PI * sqrtS);
}

/// Two body phsp factor.
/// From PDG2014 Eqn.47-2
/// \param sqrtS invariant mass of particles A and B
/// \param ma Mass of particle A
/// \param mb Mass of particle B
inline std::complex<double> phspFactorAC(double sqrtS, double ma, double mb) {
  // The PDG document states that these formulas are only in case of ma = mb
  // and the implementation is also incorrect. Use this function at own risk!
  double s = sqrtS * sqrtS;
  std::complex<double> i(0, 1);

  // == Two types of analytic continuation
  // 1) Complex sqrt
  //  std::complex<double> rhoOld;
  //  rhoOld = sqrt(std::complex<double>(qSqValue(sqrtS,ma,mb))) /
  //(8*M_PI*sqrtS); //PDG definition
  //  rhoOld = sqrt(std::complex<double>(qSqValue(sqrtS,ma,mb))) /
  //(0.5*sqrtS); //BaBar definition
  //  return rhoOld; //complex sqrt

  // 2) Correct analytic continuation
  // proper analytic continuation (PDG 2014 - Resonances (47.2.2))
  // I'm not sure of this is correct for the case of two different masses ma
  // and mb. Furthermore we divide by the factor 16*Pi*Sqrt[s]). This is
  // more or less arbitrary and not mentioned in the reference, but it leads
  // to a good agreement between both approaches.
  std::complex<double> rho;
  double q = std::sqrt(std::fabs(qSquared(s, ma, mb) * 4 / s));
  if ((ma + mb) * (ma + mb) < s) { // above threshold
    rho = (-q / M_PI * log(std::fabs((1 + q) / (1 - q))) + i * q) /
          (i * 16.0 * M_PI * sqrtS);
  } else if (0 < s && s <= (ma + mb) * (ma + mb)) { // below threshold
    rho = (-2.0 * q / M_PI * atan(1 / q)) / (i * 16.0 * M_PI * sqrtS);
  } else if (s <= 0) { // below 0
    rho = -q / M_PI * std::log(std::fabs((1 + q) / (1 - q)));
  } else
    throw std::runtime_error("phspFactor() | PhspFactor not "
                             "defined at sqrtS = " +
                             std::to_string((long double)sqrtS));

#ifndef NDEBUG
  if (rho.real() != rho.real() || rho.imag() != rho.imag()) {
    std::stringstream ss;
    ss << "phspFactor() | Result invalid (NaN)! sqrtS=" << sqrtS
       << ", ma=" << ma << ", mb=" << mb;
    throw std::runtime_error(ss.str());
  }
#endif

  return rho; // correct analytical continuation
}

/// Calculate Break-up momentum.
/// At energy \p sqrtS for particles with masses \p ma and \p mb. From
/// PDG2014 Eq.46-20a. Below threshold the function is analytically continued.
/// \param sqrtS center-of-mass energy
/// \param ma mass particle A
/// \param mb mass particle B
inline std::complex<double> qValueAC(double sqrtS, double ma, double mb) {
  return phspFactorAC(sqrtS, ma, mb) * 8.0 * M_PI * sqrtS;
}

/// Defines interface for form factors
/// It should be noted that when exchanging various form factor implementations
/// in the code, no correctness of the mathematical description is guaranteed.
class FormFactor {
public:
  virtual double operator()(double QSquared, unsigned int L,
                            double MesonRadius) const = 0;
  virtual std::string getName() const = 0;
};

class NoFormFactor : public FormFactor {
public:
  double operator()(double QSquared, unsigned int L, double MesonRadius) const {
    return 1.0;
  }
  std::string getName() const { return "NoFormFactor"; }
};

/// Blatt-Weisskopf form factors with normalization F(x=mR) = 1.
/// Reference: S.U.Chung Annalen der Physik 4(1995) 404-430
/// z = q / (interaction range). For the interaction range we assume
/// 1/mesonRadius
class BlattWeisskopfFormFactor : public FormFactor {
  double operator()(double qSquared, unsigned int L, double MesonRadius) const {
    if (MesonRadius == 0.) {
      return 1.0;
    }

    double z = qSquared * MesonRadius * MesonRadius;
    // Events below threshold. What should we do if event is below threshold?
    // Shouldn't really influence the result because resonances at threshold
    // don't have spin(?).
    // Ref. for Blatt-Weisskopf: Phys. Rev. D 48, 1225
    z = std::abs(z);

    double ff(1.0);
    if (L == 1) {
      ff = std::sqrt(2 * z / (z + 1));
    } else if (L == 2) {
      ff = std::sqrt(13 * z * z / ((z - 3) * (z - 3) + 9 * z));
    } else if (L == 3) {
      ff = std::sqrt(277 * z * z * z /
                     (z * (z - 15) * (z - 15) + 9 * (2 * z - 5) * (2 * z - 5)));
    } else if (L == 4) {
      ff = std::sqrt(12746 * z * z * z * z /
                     ((z * z - 45 * z + 105) * (z * z - 45 * z + 105) +
                      25 * z * (2 * z - 21) * (2 * z - 21)));
    }

    return ff;
  }
  std::string getName() const { return "BlattWeisskopf"; }
};

/// Form factor for a0(980) used by Crystal Barrel (Phys.Rev.D78-074023)
class CrystalBarrelFormFactor : public FormFactor {
  double operator()(double qSquared, unsigned int L, double MesonRadius) const {
    if (MesonRadius == 0.) {
      return 1.0;
    }

    double ff(1.0);
    if (L == 0) {
      double alpha = MesonRadius * MesonRadius / 6;
      ff = std::exp(-alpha * qSquared);
    }

    return ff;
  }
  std::string getName() const { return "CrystalBarrel"; }
};

std::shared_ptr<ComPWA::FunctionTree::TreeNode> createFunctionTree(
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquaredDaughter1,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquaredDaughter2,
    std::shared_ptr<ComPWA::FunctionTree::FitParameter> MesonRadius,
    unsigned int L, std::shared_ptr<FormFactor> FormFactorFunctor,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquared);

class FormFactorStrategy : public ComPWA::FunctionTree::Strategy {
public:
  FormFactorStrategy(std::shared_ptr<FormFactor> FF)
      : ComPWA::FunctionTree::Strategy(ComPWA::FunctionTree::ParType::MDOUBLE,
                                       "FormFactor"),
        FormFactorFunctor(FF) {}

  virtual void execute(ComPWA::FunctionTree::ParameterList &paras,
                       std::shared_ptr<ComPWA::FunctionTree::Parameter> &out);

private:
  std::shared_ptr<FormFactor> FormFactorFunctor;
};

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
