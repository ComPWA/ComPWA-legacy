// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
#include <iterator>
#include <numeric>
#include <vector>

#include "Coupling.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "RelativisticBreitWigner.hpp"

#include <boost/property_tree/ptree.hpp>

namespace ComPWA {
namespace Physics {
namespace Dynamics {

RelativisticBreitWigner::RelativisticBreitWigner(
    std::string name, std::pair<std::string, std::string> daughters,
    std::shared_ptr<ComPWA::PartList> partL)
    : AbstractDynamicalFunction(name) {

  LOG(TRACE) << "RelativisticBreitWigner::Factory() | Construction of " << name
             << ".";

  auto partProp = partL->find(name)->second;
  Mass = std::make_shared<FitParameter>(partProp.GetMassPar());

  auto decayTr = partProp.GetDecayInfo();
  if (partProp.GetDecayType() != "relativisticBreitWigner")
    throw std::runtime_error(
        "RelativisticBreitWigner::Factory() | Decay type does not match! ");

  auto spin = partProp.GetSpinQuantumNumber("Spin");
  J = spin;
  // in default, using spin J as Orbital Angular Momentum
  // update by calling SetOrbitalAngularMomentum() before any further process
  // after RelBW is created by calling of constructor
  L = spin;

  FFType = FormFactorType(decayTr.get<int>("FormFactor.<xmlattr>.Type"));

  // Read parameters from tree. Currently parameters of type 'Width' and
  // 'MesonRadius' are required.
  for (const auto &v : decayTr.get_child("")) {
    if (v.first != "Parameter")
      continue;
    std::string type = v.second.get<std::string>("<xmlattr>.Type");
    if (type == "Width") {
      SetWidthParameter(std::make_shared<FitParameter>(v.second));
    } else if (type == "MesonRadius") {
      SetMesonRadiusParameter(std::make_shared<FitParameter>(v.second));
    } else {
      throw std::runtime_error(
          "RelativisticBreitWigner::Factory() | Parameter of type " + type +
          " is unknown.");
    }
  }

  DaughterMasses =
      std::make_pair(partL->find(daughters.first)->second.GetMass(),
                     partL->find(daughters.second)->second.GetMass());

  DaughterNames = daughters;

  LOG(TRACE)
      << "RelativisticBreitWigner::Factory() | Construction of the decay "
      << partProp.name() << " -> " << daughters.first << " + "
      << daughters.second;
}

RelativisticBreitWigner::~RelativisticBreitWigner() {}

std::complex<double> RelativisticBreitWigner::evaluate(const DataPoint &point,
                                                       unsigned int pos) const {
  std::complex<double> result = dynamicalFunction(
      point.KinematicVariableList[pos], Mass->value(), DaughterMasses.first,
      DaughterMasses.second, Width->value(), (double)L, MesonRadius->value(),
      FFType);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return result;
}

std::complex<double> RelativisticBreitWigner::dynamicalFunction(
    double mSq, double mR, double ma, double mb, double width, unsigned int L,
    double mesonRadius, ComPWA::Physics::Dynamics::FormFactorType ffType) {

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
  // Barrier term (PDG 2014 Eq. 47.23)
  // \f[
  //     barrierTermSq = \left( \frac{q(s)}{q(s_R)} \right)^{2L+1} \times
  //                     \left( \frac{F(s)}{F(s_R)} \right)^{2}
  // \f]
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

std::shared_ptr<ComPWA::FunctionTree>
RelativisticBreitWigner::createFunctionTree(const ParameterList &DataSample,
                                            unsigned int pos,
                                            const std::string &suffix) const {

  // size_t sampleSize = DataSample.mDoubleValue(pos)->values().size();
  size_t sampleSize = DataSample.mDoubleValue(0)->values().size();

  auto tr = std::make_shared<FunctionTree>(
      "RelBreitWigner" + suffix, MComplex("", sampleSize),
      std::make_shared<BreitWignerStrategy>());

  tr->createLeaf("Mass", Mass, "RelBreitWigner" + suffix);
  tr->createLeaf("Width", Width, "RelBreitWigner" + suffix);
  tr->createLeaf("OrbitalAngularMomentum", (double)L,
                 "RelBreitWigner" + suffix);
  tr->createLeaf("MesonRadius", MesonRadius, "RelBreitWigner" + suffix);
  tr->createLeaf("FormFactorType", FFType, "RelBreitWigner" + suffix);
  tr->createLeaf("MassA", DaughterMasses.first, "RelBreitWigner" + suffix);
  tr->createLeaf("MassB", DaughterMasses.second, "RelBreitWigner" + suffix);
  /*tr->createLeaf("Data_mSq[" + std::to_string(pos) + "]",
                  DataSample.mDoubleValue(pos), "RelBreitWigner" + suffix);*/
  tr->createLeaf("Data_mSq[" + std::to_string(pos) + "]",
                 DataSample.mDoubleValue(pos), "RelBreitWigner" + suffix);

  return tr;
};

void BreitWignerStrategy::execute(ParameterList &paras,
                                  std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter(
        "BreitWignerStrat::execute() | Parameter type mismatch!");

#ifndef NDEBUG
  // Check parameter type
  if (checkType != out->type())
    throw(WrongParType("BreitWignerStrat::execute() | "
                       "Output parameter is of type " +
                       std::string(ParNames[out->type()]) +
                       " and conflicts with expected type " +
                       std::string(ParNames[checkType])));

  // How many parameters do we expect?
  size_t check_nInt = 0;
  size_t nInt = paras.intValues().size();
  size_t check_nDouble = 7;
  size_t nDouble = paras.doubleValues().size();
  nDouble += paras.doubleParameters().size();
  size_t check_nComplex = 0;
  size_t nComplex = paras.complexValues().size();
  size_t check_nMInteger = 0;
  size_t nMInteger = paras.mIntValues().size();
  size_t check_nMDouble = 1;
  size_t nMDouble = paras.mDoubleValues().size();
  size_t check_nMComplex = 0;
  size_t nMComplex = paras.mComplexValues().size();

  // Check size of parameter list
  if (nInt != check_nInt)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of IntParameters does not match: " +
                       std::to_string(nInt) + " given but " +
                       std::to_string(check_nInt) + " expected."));
  if (nDouble != check_nDouble)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of FitParameters does not match: " +
                       std::to_string(nDouble) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (nComplex != check_nComplex)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(nComplex) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (nMInteger != check_nMInteger)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of MultiInt does not match: " +
                       std::to_string(nMInteger) + " given but " +
                       std::to_string(check_nMInteger) + " expected."));
  if (nMDouble != check_nMDouble)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of MultiDoubles does not match: " +
                       std::to_string(nMDouble) + " given but " +
                       std::to_string(check_nMDouble) + " expected."));
  if (nMComplex != check_nMComplex)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of MultiComplexes does not match: " +
                       std::to_string(nMComplex) + " given but " +
                       std::to_string(check_nMComplex) + " expected."));
#endif

  size_t n = paras.mDoubleValue(0)->values().size();
  if (!out)
    out = MComplex("", n);
  auto par =
      std::static_pointer_cast<Value<std::vector<std::complex<double>>>>(out);
  auto &results = par->values(); // reference
  if (results.size() != n) {
    results.resize(n);
  }
  // Get parameters from ParameterList:
  // We use the same order of the parameters as was used during tree
  // construction.
  double m0 = paras.doubleParameter(0)->value();
  double Gamma0 = paras.doubleParameter(1)->value();
  double MesonRadius = paras.doubleParameter(2)->value();
  unsigned int orbitL = paras.doubleValue(0)->value();
  FormFactorType ffType = FormFactorType(paras.doubleValue(1)->value());
  double ma = paras.doubleValue(2)->value();
  double mb = paras.doubleValue(3)->value();

  // calc function for each point
  for (unsigned int ele = 0; ele < n; ele++) {
    try {
      results.at(ele) = RelativisticBreitWigner::dynamicalFunction(
          paras.mDoubleValue(0)->values().at(ele), m0, ma, mb, Gamma0, orbitL,
          MesonRadius, ffType);
    } catch (std::exception &ex) {
      LOG(ERROR) << "BreitWignerStrategy::execute() | " << ex.what();
      throw(std::runtime_error("BreitWignerStrategy::execute() | "
                               "Evaluation of dynamic function failed!"));
    }
  }
}

void RelativisticBreitWigner::addUniqueParametersTo(ParameterList &list) {
  // We check of for each parameter if a parameter of the same name exists in
  // list. If so we check if both are equal and set the local parameter to the
  // parameter from the list. In this way we connect parameters that occur on
  // different positions in the amplitude.
  Mass = list.addUniqueParameter(Mass);
  Width = list.addUniqueParameter(Width);
  MesonRadius = list.addUniqueParameter(MesonRadius);
}

void RelativisticBreitWigner::addFitParametersTo(std::vector<double> &FitParameters) {
  FitParameters.push_back(Mass->value());
  FitParameters.push_back(Width->value());
  FitParameters.push_back(MesonRadius->value());
}

void RelativisticBreitWigner::updateParametersFrom(const ParameterList &list) {
  // Try to update Mass
  std::shared_ptr<FitParameter> mass;
  try {
    mass = FindParameter(Mass->name(), list);
  } catch (std::exception &ex) {
  }
  if (mass)
    Mass->updateParameter(mass);

  // Try to update mesonRadius
  std::shared_ptr<FitParameter> rad;
  try {
    rad = FindParameter(MesonRadius->name(), list);
  } catch (std::exception &ex) {
  }
  if (rad)
    MesonRadius->updateParameter(rad);

  // Try to update width
  std::shared_ptr<FitParameter> width;
  try {
    width = FindParameter(Width->name(), list);
  } catch (std::exception &ex) {
  }
  if (width)
    Width->updateParameter(width);

  return;
}

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA
