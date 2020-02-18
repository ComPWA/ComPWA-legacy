// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_DYNAMICS_FLATTE_HPP_
#define COMPWA_PHYSICS_DYNAMICS_FLATTE_HPP_

#include "Coupling.hpp"
#include "FormFactor.hpp"
#include "RelativisticBreitWigner.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

namespace Flatte {
struct InputInfo : ComPWA::Physics::Dynamics::InputInfo {
  /// Coupling to signal channel
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> G;
  /// Coupling parameters and final state masses for multiple hidden channels
  std::vector<Coupling> HiddenCouplings;
};

/// Helper function to calculate the coupling terms for the Flatte formular.
inline std::complex<double>
flatteCouplingTerm(double sqrtS, double mR, double coupling, double massA,
                   double massB, unsigned int J, double mesonRadius,
                   std::shared_ptr<FormFactor> FormFactorFunctor) {
  auto qR = std::abs(qValueAC(mR, massA, massB));
  auto phspR = phspFactorAC(sqrtS, massA, massB);
  auto ffR = FormFactorFunctor->operator()(qR *qR, J, mesonRadius);
  auto qS = std::abs(qValueAC(sqrtS, massA, massB));
  auto barrierA = FormFactorFunctor->operator()(qS *qS, J, mesonRadius) / ffR;

  // Calculate normalized vertex functions vtxA(s_R)
  std::complex<double> vtxA(1, 0); // spin==0
  if (J > 0 || FormFactorFunctor->getName() == "CrystalBarrel") {
    vtxA = ffR * std::pow(qR, J);
  }
  auto width = couplingToWidth(mR, coupling, vtxA, phspR);
  // Including the factor qTermA, as suggested by PDG 2014, Chapter 47.2,
  // leads to an amplitude that doesn't converge.
  //  qTermA = qValue(sqrtS,massA1,massA2) / qValue(mR,massA1,massA2);
  //  termA = gammaA * barrierA * barrierA * std::pow(qTermA, (double)2 * J +
  //  1);

  return (width * barrierA * barrierA);
}

/** Dynamical function for two coupled channel approach
 *
 * @param mSq center-of-mass energy^2 (=s)
 * @param mR mass of resonances
 * @param gA coupling constant for signal channel
 * @param termA Coupling term to signal channel
 * @param termB Coupling term to second channel
 * @param termC Coupling term to third channel (optional)
 * @return
 */
inline std::complex<double>
dynamicalFunction(double mSq, double mR, double gA, std::complex<double> termA,
                  std::complex<double> termB,
                  std::complex<double> termC = std::complex<double>(0, 0)) {
  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  std::complex<double> denom = std::complex<double>(mR * mR - mSq, 0) +
                               (-1.0) * i * sqrtS * (termA + termB + termC);

  std::complex<double> result = std::complex<double>(gA, 0) / denom;

#ifndef NDEBUG
  if (std::isnan(result.real()) || std::isnan(result.imag())) {
    std::cout << "AmpFlatteRes::dynamicalFunction() | " << mR << " " << mSq
              << " " << termA << " " << termB << " " << termC << std::endl;
    return 0;
  }
#endif

  return result;
}

/** Dynamical function for two coupled channel approach
 *
 * @param mSq center-of-mass energy^2 (=s)
 * @param mR mass of resonances
 * @param massA1 mass of first particle of signal channel
 * @param massA2 mass of second particle of signal channel
 * @param gA coupling constant for signal channel
 * @param massB1 mass of first particle of second channel
 * @param massB2 mass of second particle of second channel
 * @param couplingB coupling constant for second channel
 * @param massC1 mass of first particle of third channel
 * @param massC2 mass of third particle of third channel
 * @param couplingC coupling constant for third channel
 * @param L Orbital angular momentum between two daughters a and b
 * @param mesonRadius 1/interaction length (needed for barrier factors)
 * @param FormFactorFunctor functor of the form factor
 * @return
 */
inline std::complex<double>
dynamicalFunction(double mSq, double mR, double massA1, double massA2,
                  double gA, double massB1, double massB2, double couplingB,
                  double massC1, double massC2, double couplingC,
                  unsigned int L, double mesonRadius,
                  std::shared_ptr<FormFactor> FormFactorFunctor) {
  double sqrtS = sqrt(mSq);

  // channel A - signal channel
  auto termA = flatteCouplingTerm(sqrtS, mR, gA, massA1, massA2, L, mesonRadius,
                                  FormFactorFunctor);
  // channel B - hidden channel
  auto termB = flatteCouplingTerm(sqrtS, mR, couplingB, massB1, massB2, L,
                                  mesonRadius, FormFactorFunctor);

  // channel C - hidden channel
  std::complex<double> termC;
  if (couplingC != 0.0) {
    termC = flatteCouplingTerm(sqrtS, mR, couplingC, massC1, massC2, L,
                               mesonRadius, FormFactorFunctor);
  }
  return dynamicalFunction(mSq, mR, gA, termA, termB, termC);
}

std::shared_ptr<ComPWA::FunctionTree::TreeNode> createFunctionTree(
    InputInfo Params,
    std::shared_ptr<ComPWA::FunctionTree::Value<std::vector<double>>>
        InvMassSquared);

} // namespace Flatte

class FlatteStrategy : public ComPWA::FunctionTree::Strategy {
public:
  FlatteStrategy(std::shared_ptr<FormFactor> FF)
      : Strategy(ComPWA::FunctionTree::ParType::MCOMPLEX, "Flatte"),
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
