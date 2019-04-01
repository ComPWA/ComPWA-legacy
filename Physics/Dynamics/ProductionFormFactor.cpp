// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
#include <iterator>
#include <numeric>
#include <vector>

#include "Coupling.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "ProductionFormFactor.hpp"

#include <boost/property_tree/ptree.hpp>

namespace ComPWA {
namespace Physics {
namespace Dynamics {

ProductionFormFactor::ProductionFormFactor(
    std::string name, std::pair<std::string, std::string> daughters,
    std::shared_ptr<ComPWA::PartList> partL)
    : AbstractDynamicalFunction(name) {

  LOG(TRACE) << "ProductionFormFactor::Factory() | Construction of " << name
             << ".";

  auto partProp = partL->find(name)->second;

  auto decayTr = partProp.GetDecayInfo();

  auto spin = partProp.GetSpinQuantumNumber("Spin");
  J = spin;
  // in default, using spin J as Orbital Angular Momentum
  // update by calling SetOrbitalAngularMomentum() before any further process
  // after the ProductionFormFactor is created by calling of constructor
  L = spin;

  FFType = FormFactorType(decayTr.get<int>("FormFactor.<xmlattr>.Type"));

  // Read parameters from tree. Currently parameters of type 'MesonRadius' 
  // is required.
  // But width is also in this tree, because it is needed by RelBW
  for (const auto &v : decayTr.get_child("")) {
    if ("Parameter" != v.first)
      continue;
    std::string type = v.second.get<std::string>("<xmlattr>.Type");
    if ("MesonRadius" == type) {
      SetMesonRadiusParameter(std::make_shared<FitParameter>(v.second));
    } else if ("Width" == type) {
      continue;
    } else {
      throw std::runtime_error(
          "ProductionFormFactor::Factory() | Parameter of type " + type +
          " is unknown.");
    }
  }

  DaughterMasses =
      std::make_pair(partL->find(daughters.first)->second.GetMass(),
                     partL->find(daughters.second)->second.GetMass());

  DaughterNames = daughters;

  LOG(TRACE)
      << "ProductionFormFactor::Factory() | Construction of the decay "
      << partProp.name() << " -> " << daughters.first << " + "
      << daughters.second;
}

ProductionFormFactor::~ProductionFormFactor() {}

std::complex<double> ProductionFormFactor::evaluate(const DataPoint &point,
                                                       unsigned int pos) const {
  std::complex<double> result = dynamicalFunction(
      point.KinematicVariableList[pos], DaughterMasses.first, DaughterMasses.second, 
      (unsigned int)L, MesonRadius->value(), FFType);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return result;
}

std::complex<double> ProductionFormFactor::dynamicalFunction(
    double mSq, double ma, double mb, unsigned int L, double mesonRadius, 
    ComPWA::Physics::Dynamics::FormFactorType ffType) {

  std::complex<double> i(0, 1);
  double sqrtS = std::sqrt(mSq);
  
  // currently call ComPWA::Physics::Dynamics::FormFactor() to calculate the
  // form factor. In furthure, we may use a FormFactor class to provide both
  // production and decay form factor and merge ProductionFormFactor class 
  // and FormFactor functions
  double ff = FormFactor(sqrtS, ma, mb, L, mesonRadius, ffType); 
  
  return std::complex<double>(ff, 0);
}

std::shared_ptr<ComPWA::FunctionTree>
ProductionFormFactor::createFunctionTree(const ParameterList &DataSample,
                                            unsigned int pos,
                                            const std::string &suffix) const {

  // size_t sampleSize = DataSample.mDoubleValue(pos)->values().size();
  size_t sampleSize = DataSample.mDoubleValue(0)->values().size();

  auto tr = std::make_shared<FunctionTree>(
      "ProductionFormFactor" + suffix, MComplex("", sampleSize),
      std::make_shared<FormFactorStrategy>());

  tr->createLeaf("OrbitalAngularMomentum", (unsigned int)L,
                 "ProductionFormFactor" + suffix);
  tr->createLeaf("MesonRadius", MesonRadius, "ProductionFormFactor" + suffix);
  tr->createLeaf("FormFactorType", FFType, "ProductionFormFactor" + suffix);
  tr->createLeaf("MassA", DaughterMasses.first, "ProductionFormFactor" + suffix);
  tr->createLeaf("MassB", DaughterMasses.second, "ProductionFormFactor" + suffix);
  tr->createLeaf("Data_mSq[" + std::to_string(pos) + "]",
                 DataSample.mDoubleValue(pos), "ProductionFormFactor" + suffix);

  return tr;
};

void FormFactorStrategy::execute(ParameterList &paras,
                                  std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter(
        "FormFactorStrat::execute() | Parameter type mismatch!");

#ifndef NDEBUG
  // Check parameter type
  if (checkType != out->type())
    throw(WrongParType("FormFactorStrat::execute() | "
                       "Output parameter is of type " +
                       std::string(ParNames[out->type()]) +
                       " and conflicts with expected type " +
                       std::string(ParNames[checkType])));

  // How many parameters do we expect?
  size_t check_nInt = 2;
  size_t nInt = paras.intValues().size();
  size_t check_nDouble = 3;
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
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of IntParameters does not match: " +
                       std::to_string(nInt) + " given but " +
                       std::to_string(check_nInt) + " expected."));
  if (nDouble != check_nDouble)
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of FitParameters does not match: " +
                       std::to_string(nDouble) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (nComplex != check_nComplex)
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(nComplex) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (nMInteger != check_nMInteger)
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of MultiInt does not match: " +
                       std::to_string(nMInteger) + " given but " +
                       std::to_string(check_nMInteger) + " expected."));
  if (nMDouble != check_nMDouble)
    throw(BadParameter("FormFactorStrat::execute() | "
                       "Number of MultiDoubles does not match: " +
                       std::to_string(nMDouble) + " given but " +
                       std::to_string(check_nMDouble) + " expected."));
  if (nMComplex != check_nMComplex)
    throw(BadParameter("FormFactorStrat::execute() | "
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
  double MesonRadius = paras.doubleParameter(0)->value();
  unsigned int orbitL = paras.intValue(0)->value();
  FormFactorType ffType = FormFactorType(paras.intValue(1)->value());
  double ma = paras.doubleValue(1)->value();
  double mb = paras.doubleValue(2)->value();

  // calc function for each point
  for (unsigned int ele = 0; ele < n; ele++) {
    try {
      results.at(ele) = ProductionFormFactor::dynamicalFunction(
          paras.mDoubleValue(0)->values().at(ele), ma, mb, orbitL,
          MesonRadius, ffType);
    } catch (std::exception &ex) {
      LOG(ERROR) << "FormFactorStrategy::execute() | " << ex.what();
      throw(std::runtime_error("FormFactorStrategy::execute() | "
                               "Evaluation of dynamic function failed!"));
    }
  }
}

void ProductionFormFactor::addUniqueParametersTo(ParameterList &list) {
  // We check of for each parameter if a parameter of the same name exists in
  // list. If so we check if both are equal and set the local parameter to the
  // parameter from the list. In this way we connect parameters that occur on
  // different positions in the amplitude.
  MesonRadius = list.addUniqueParameter(MesonRadius);
}

void ProductionFormFactor::addFitParametersTo(std::vector<double> &FitParameters) {
  FitParameters.push_back(MesonRadius->value());
}

void ProductionFormFactor::updateParametersFrom(const ParameterList &list) {

  // Try to update mesonRadius
  std::shared_ptr<FitParameter> rad;
  try {
    rad = FindParameter(MesonRadius->name(), list);
  } catch (std::exception &ex) {
  }
  if (rad)
    MesonRadius->updateParameter(rad);

  return;
}

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA
