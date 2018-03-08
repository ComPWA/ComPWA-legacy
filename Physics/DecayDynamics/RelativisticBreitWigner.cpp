// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
#include <numeric>
#include <iterator>
#include <vector>

#include <boost/property_tree/ptree.hpp>

#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/DecayDynamics/RelativisticBreitWigner.hpp"
#include "Physics/DecayDynamics/Coupling.hpp"

using namespace ComPWA::Physics::DecayDynamics;

RelativisticBreitWigner::RelativisticBreitWigner(
   std::string name, std::pair<std::string,std::string> daughters,
               std::shared_ptr<ComPWA::PartList> partL) {

  LOG(trace) << "RelativisticBreitWigner::Factory() | Construction of " << name
             << ".";
  setName(name);
  auto partProp = partL->find(name)->second;
  SetMassParameter(
      std::make_shared<FitParameter>(partProp.GetMassPar()));

  auto decayTr = partProp.GetDecayInfo();
  if (partProp.GetDecayType() != "relativisticBreitWigner")
    throw std::runtime_error(
        "RelativisticBreitWigner::Factory() | Decay type does not match! ");

  auto spin = partProp.GetSpinQuantumNumber("Spin");
  SetSpin(spin);

  auto ffType = formFactorType(decayTr.get<int>("FormFactor.<xmlattr>.Type"));
  SetFormFactorType(ffType);

  // Read parameters from tree. Currently parameters of type 'Width' and
  // 'MesonRadius' are required.
  for (const auto &v : decayTr.get_child("")) {
    if (v.first != "Parameter")
      continue;
    std::string type = v.second.get<std::string>("<xmlattr>.Type");
    if (type == "Width") {
      SetWidthParameter(std::make_shared<FitParameter>(v.second));
    } else if (type == "MesonRadius") {
      SetMesonRadiusParameter(
          std::make_shared<FitParameter>(v.second));
    } else {
      throw std::runtime_error(
          "RelativisticBreitWigner::Factory() | Parameter of type " + type +
          " is unknown.");
    }
  }

  std::pair<double, double> daughterMasses(
      partL->find(daughters.first)->second.GetMass(),
      partL->find(daughters.second)->second.GetMass());
  SetDecayMasses(daughterMasses);
  SetDecayNames(daughters);

  LOG(trace)
      << "RelativisticBreitWigner::Factory() | Construction of the decay "
      << partProp.name() << " -> " << daughters.first << " + "
      << daughters.second;
}

std::complex<double> RelativisticBreitWigner::evaluate(const DataPoint &point,
                                                       int pos) const {
  std::complex<double> result =
      dynamicalFunction(point.value(pos),Mass->value(), DaughterMasses.first,
                        DaughterMasses.second,Width->value(), (double)J,
                        MesonRadius->value(), FormFactorType);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return result;
}

bool RelativisticBreitWigner::isModified() const {
  if (AbstractDynamicalFunction::isModified())
    return true;
  if (Width->value() != CurrentWidth ||
      MesonRadius->value() != CurrentMesonRadius) {
    setModified();
    const_cast<double &>(CurrentWidth) =Width->value();
    const_cast<double &>(CurrentMesonRadius) = MesonRadius->value();
    return true;
  }
  return false;
}

std::complex<double> RelativisticBreitWigner::dynamicalFunction(
    double mSq, double mR, double ma, double mb, double width, unsigned int J,
    double mesonRadius, formFactorType ffType) {

  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  // Phase space factors at sqrt(s) and at the resonance position
  auto phspFactorSqrtS = phspFactor(sqrtS, ma, mb);
  auto phspFactormR = phspFactor(mR, ma, mb);

  // Check if we have an event which is exactly at the phase space boundary
  if (phspFactorSqrtS == std::complex<double>(0, 0))
    return std::complex<double>(0, 0);

  std::complex<double> qRatio =
      std::pow((phspFactorSqrtS / phspFactormR) * mR / sqrtS, (2 * J + 1));
  double ffR = FormFactor(mR, ma, mb, J, mesonRadius, ffType);
  // Barrier factor ( F(sqrt(s)) / F(mR) )
  double barrier = FormFactor(sqrtS, ma, mb, J, mesonRadius, ffType) / ffR;
  std::complex<double> damping = barrier*qRatio;

  // Calculate normalized vertex function gammaA(s_R) (see PDG2014, Chapter 47.2)
  std::complex<double> gammaA(1, 0); // spin==0
  if (J > 0) {
    std::complex<double> qR = std::pow(qValue(mR, ma, mb), J);
    gammaA = ffR * qR;
  }
  
  // Coupling to the final state (ma, mb)
  std::complex<double> g_final = widthToCoupling(mR, width, gammaA, phspFactorSqrtS);

  std::complex<double> denom(mR * mR - mSq, 0);
  denom += (-1.0) * i * sqrtS * (width * damping);

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
RelativisticBreitWigner::tree(const ParameterList &sample, int pos,
                                 std::string suffix) {

  size_t sampleSize = sample.mDoubleValue(pos)->values().size();

  auto tr = std::make_shared<FunctionTree>(
      "RelBreitWigner" + suffix, MComplex("", sampleSize),
      std::make_shared<BreitWignerStrategy>());

  tr->createLeaf("Mass",Mass, "RelBreitWigner" + suffix);
  tr->createLeaf("Width",Width, "RelBreitWigner" + suffix);
  tr->createLeaf("Spin", (double)J, "RelBreitWigner" + suffix);
  tr->createLeaf("MesonRadius", MesonRadius, "RelBreitWigner" + suffix);
  tr->createLeaf("FormFactorType", FormFactorType, "RelBreitWigner" + suffix);
  tr->createLeaf("MassA", DaughterMasses.first, "RelBreitWigner" + suffix);
  tr->createLeaf("MassB", DaughterMasses.second, "RelBreitWigner" + suffix);
  tr->createLeaf("Data_mSq[" + std::to_string(pos) + "]",
                 sample.mDoubleValue(pos), "RelBreitWigner" + suffix);

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

  // Get parameters from ParameterList:
  // We use the same order of the parameters as was used during tree
  // construction.
  double m0 = paras.doubleParameter(0)->value();
  double Gamma0 = paras.doubleParameter(1)->value();
  double d = paras.doubleParameter(2)->value();
  unsigned int spin = paras.doubleValue(0)->value();
  formFactorType ffType = formFactorType(paras.doubleValue(1)->value());
  double ma = paras.doubleValue(2)->value();
  double mb = paras.doubleValue(3)->value();

  // calc function for each point
  for (unsigned int ele = 0; ele < n; ele++) {
    try {
      results.at(ele) = RelativisticBreitWigner::dynamicalFunction(
          paras.mDoubleValue(0)->values().at(ele), m0, ma, mb, Gamma0, spin, d,
          ffType);
    } catch (std::exception &ex) {
      LOG(error) << "BreitWignerStrategy::execute() | " << ex.what();
      throw(std::runtime_error("BreitWignerStrategy::execute() | "
                               "Evaluation of dynamic function failed!"));
    }
  }
}

void RelativisticBreitWigner::parameters(ParameterList &list) {
  AbstractDynamicalFunction::parameters(list);

  // We check of for each parameter if a parameter of the same name exists in
  // list. If so we check if both are equal and set the local parameter to the
  // parameter from the list. In this way we connect parameters that occur on
  // different positions in the amplitude.
 Width = list.addUniqueParameter(Width);
  MesonRadius = list.addUniqueParameter(MesonRadius);
}

void RelativisticBreitWigner::updateParameters(const ParameterList &list) {

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
