// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_DYNAMICS_FORMFACTOR_HPP_
#define COMPWA_PHYSICS_DYNAMICS_FORMFACTOR_HPP_

#include <complex>
#include <exception>

namespace ComPWA {
namespace Physics {
namespace Dynamics {

/// Calculate Break-up momentum squared.
/// At energy \p sqrtS for particles with masses \p ma and \p mb. From
/// PDG2014 Eq.46-20a. Below threshold the function is analytically continued.
/// \param sqrtS center-of-mass energy
/// \param ma mass particle A
/// \param mb mass particle B
inline double qSqValue(double sqrtS, double ma, double mb) {
  double mapb = ma + mb;
  double mamb = ma - mb;
  double xSq = sqrtS * sqrtS;
  double t1 = xSq - mapb * mapb;
  double t2 = xSq - mamb * mamb;
  return (t1 * t2 / (4 * xSq));
}

/// Two body phsp factor.
/// From PDG2014 Eqn.47-2
/// \param sqrtS invariant mass of particles A and B
/// \param ma Mass of particle A
/// \param mb Mass of particle B
inline std::complex<double> phspFactor(double sqrtS, double ma, double mb) {
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
  double q = std::sqrt(std::fabs(qSqValue(sqrtS, ma, mb) * 4 / s));
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
inline std::complex<double> qValue(double sqrtS, double ma, double mb) {
  return phspFactor(sqrtS, ma, mb) * 8.0 * M_PI * sqrtS;
}

static const char *formFactorTypeString[] = {"noFormFactor", "BlattWeisskopf",
                                             "CrystalBarrel"};

enum FormFactorType { noFormFactor = 0, BlattWeisskopf = 1, CrystalBarrel = 2 };

/// Calculate form factor from the (complex) break-up momentum \p qValue.
inline double FormFactor(std::complex<double> qValue, unsigned int orbitL,
                         double mesonRadius, FormFactorType type) {
  if (mesonRadius == 0)
    return 1.0; // disable form factors
  if (type == FormFactorType::noFormFactor)
    return 1.0; // disable form factors
  if (type == FormFactorType::BlattWeisskopf && orbitL == 0) {
    return 1.0;
  }

  double ff = 0.0;

  if (type == FormFactorType::CrystalBarrel) {
    // Form factor for a0(980) used by Crystal Barrel (Phys.Rev.D78-074023)
    if (orbitL == 0) {
      double qSq = std::norm(qValue);
      double alpha = mesonRadius * mesonRadius / 6;
      ff = std::exp(-alpha * qSq);
    } else
      throw std::runtime_error("FormFactor() | Form factors of type " +
                               std::string(formFactorTypeString[type]) +
                               " are implemented for spin 0 only!");
  } else if (type == FormFactorType::BlattWeisskopf) {
    // Blatt-Weisskopt form factors with normalization F(x=mR) = 1.
    // Reference: S.U.Chung Annalen der Physik 4(1995) 404-430
    // z = q / (interaction range). For the interaction range we assume
    // 1/mesonRadius
    if (orbitL == 0)
      return 1.0;

    double qSq = std::norm(qValue);
    double z = qSq * mesonRadius * mesonRadius;
    // Events below threshold. What should we do if event is below threshold?
    // Shouldn't really influence the result because resonances at threshold
    // don't have spin(?).
    // Ref. for Blatt-Weisskopf: Phys. Rev. D 48, 1225
    z = std::fabs(z);

    if (orbitL == 1) {
      ff = std::sqrt(2 * z / (z + 1));
    } else if (orbitL == 2) {
      ff = std::sqrt(13 * z * z / ((z - 3) * (z - 3) + 9 * z));
    } else if (orbitL == 3) {
      ff = std::sqrt(277 * z * z * z /
                     (z * (z - 15) * (z - 15) + 9 * (2 * z - 5) * (2 * z - 5)));
    } else if (orbitL == 4) {
      ff = std::sqrt(12746 * z * z * z * z /
                     ((z * z - 45 * z + 105) * (z * z - 45 * z + 105) +
                      25 * z * (2 * z - 21) * (2 * z - 21)));
    } else
      throw std::runtime_error("FormFactor() | Form factors of type " +
                               std::string(formFactorTypeString[type]) +
                               " are implemented for spins up to 4!");
  } else {
    throw std::runtime_error("FormFactor() | Form factor type " +
                             std::to_string((long long int)type) +
                             " not specified!");
  }

  return ff;
}

/// Calculate form factor from sqrt(s) and masses of the final state particles.
inline double FormFactor(double sqrtS, double ma, double mb, unsigned int orbitL,
                         double mesonRadius, FormFactorType type) {
  if (type == FormFactorType::noFormFactor) {
    return 1.0;
  }
  if (type == FormFactorType::BlattWeisskopf && orbitL == 0) {
    return 1.0;
  }

  std::complex<double> qV = qValue(sqrtS, ma, mb);

  return FormFactor(qV, orbitL, mesonRadius, type);
}

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
