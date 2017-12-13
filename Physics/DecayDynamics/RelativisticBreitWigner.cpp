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

std::shared_ptr<AbstractDynamicalFunction>
RelativisticBreitWigner::Factory(std::shared_ptr<PartList> partL,
                                 const boost::property_tree::ptree &pt) {
  auto obj = std::make_shared<RelativisticBreitWigner>();

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  LOG(trace) << "RelativisticBreitWigner::Factory() | Construction of " << name
             << ".";
  obj->setName(name);
  auto partProp = partL->find(name)->second;
  obj->SetMassParameter(
      std::make_shared<DoubleParameter>(partProp.GetMassPar()));

  auto decayTr = partProp.GetDecayInfo();
  if (partProp.GetDecayType() != "relativisticBreitWigner")
    throw std::runtime_error(
        "RelativisticBreitWigner::Factory() | Decay type does not match! ");

  auto spin = partProp.GetSpinQuantumNumber("Spin");
  obj->SetSpin(spin);

  auto ffType = formFactorType(decayTr.get<int>("FormFactor.<xmlattr>.Type"));
  obj->SetFormFactorType(ffType);

  // Read parameters from tree. Currently parameters of type 'Width' and
  // 'MesonRadius' are required.
  for (const auto &v : decayTr.get_child("")) {
    if (v.first != "Parameter")
      continue;
    std::string type = v.second.get<std::string>("<xmlattr>.Type");
    if (type == "Width") {
      auto width = DoubleParameter();
      width.load(v.second);
      obj->SetWidthParameter(std::make_shared<DoubleParameter>(width));
    } else if (type == "MesonRadius") {
      auto mesonRadius = DoubleParameter();
      mesonRadius.load(v.second);
      obj->SetMesonRadiusParameter(
          std::make_shared<DoubleParameter>(mesonRadius));
    } else {
      throw std::runtime_error(
          "RelativisticBreitWigner::Factory() | Parameter of type " + type +
          " is unknown.");
    }
  }

  // Get masses of decay products
  auto decayProducts = pt.get_child("DecayProducts");
  if (decayProducts.size() != 2)
    throw boost::property_tree::ptree_error(
        "RelativisticBreitWigner::Factory() | Expect exactly two decay "
        "products (" +
        std::to_string(decayProducts.size()) + " given)!");

  auto firstItr = decayProducts.begin();
  // auto secondItr = decayProducts.begin()+1; //compile error, no idea for why
  auto secondItr = decayProducts.begin();
  secondItr++;

  std::pair<std::string, std::string> daughterNames(
      firstItr->second.get<std::string>("<xmlattr>.Name"),
      secondItr->second.get<std::string>("<xmlattr>.Name"));
  std::pair<double, double> daughterMasses(
      partL->find(daughterNames.first)->second.GetMass(),
      partL->find(daughterNames.second)->second.GetMass());
  obj->SetDecayMasses(daughterMasses);
  obj->SetDecayNames(daughterNames);

  LOG(trace)
      << "RelativisticBreitWigner::Factory() | Construction of the decay "
      << partProp.name() << " -> " << daughterNames.first << " + "
      << daughterNames.second;

  return std::static_pointer_cast<AbstractDynamicalFunction>(obj);
}

std::complex<double> RelativisticBreitWigner::evaluate(const DataPoint &point,
                                                       int pos) const {
  std::complex<double> result =
      dynamicalFunction(point.value(pos), _mass->value(), _daughterMasses.first,
                        _daughterMasses.second, _width->value(), (double)_spin,
                        _mesonRadius->value(), _ffType);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return result;
}

bool RelativisticBreitWigner::CheckModified() const {
  if (AbstractDynamicalFunction::CheckModified())
    return true;
  if (_width->value() != _current_width ||
      _mesonRadius->value() != _current_mesonRadius) {
    SetModified();
    const_cast<double &>(_current_width) = _width->value();
    const_cast<double &>(_current_mesonRadius) = _mesonRadius->value();
    return true;
  }
  return false;
}

std::complex<double> RelativisticBreitWigner::dynamicalFunction(
    double mSq, double mR, double ma, double mb, double width, unsigned int J,
    double mesonRadius, formFactorType ffType) {

  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  auto phspFactorSqrtS = phspFactor(sqrtS, ma, mb);
  auto phspFactormR = phspFactor(mR, ma, mb);

  // Check if we have an event which is exactly at the phase space boundary
  if (phspFactorSqrtS == std::complex<double>(0, 0))
    return std::complex<double>(0, 0);

  std::complex<double> qTerm =
      std::pow((phspFactorSqrtS / phspFactormR) * mR / sqrtS, (2 * J + 1));
  double barrier = FormFactor(sqrtS, ma, mb, J, mesonRadius, ffType) /
                   FormFactor(mR, ma, mb, J, mesonRadius, ffType);

  // Calculate coupling constant to final state
  std::complex<double> g_final =
      widthToCoupling(mSq, mR, width, ma, mb, J, mesonRadius, ffType);

  // Coupling constant from production reaction. In case of a particle decay
  // the production coupling doesn't depend in energy since the CM energy
  // is in the (RC) system fixed to the mass of the decaying particle
  double g_production = 1;

  std::complex<double> denom = std::complex<double>(mR * mR - mSq, 0) +
                               (-1.0) * i * sqrtS * (width * qTerm * barrier);

  std::complex<double> result = g_final * g_production / denom;

  assert(
      (!std::isnan(result.real()) || !std::isinf(result.real())) &&
      "RelativisticBreitWigner::dynamicalFunction() | Result is NaN or Inf!");
  assert(
      (!std::isnan(result.imag()) || !std::isinf(result.imag())) &&
      "RelativisticBreitWigner::dynamicalFunction() | Result is NaN or Inf!");

  return result;
}

std::shared_ptr<ComPWA::FunctionTree>
RelativisticBreitWigner::GetTree(const ParameterList &sample, int pos,
                                 std::string suffix) {

  size_t sampleSize = sample.mDoubleValue(pos)->values().size();

  auto tr = std::make_shared<FunctionTree>(
      "RelBreitWigner" + suffix, MComplex("", sampleSize),
      std::make_shared<BreitWignerStrategy>());

  tr->createLeaf("Mass", _mass, "RelBreitWigner" + suffix);
  tr->createLeaf("Width", _width, "RelBreitWigner" + suffix);
  tr->createLeaf("Spin", (double)_spin, "RelBreitWigner" + suffix);
  tr->createLeaf("MesonRadius", _mesonRadius, "RelBreitWigner" + suffix);
  tr->createLeaf("FormFactorType", _ffType, "RelBreitWigner" + suffix);
  tr->createLeaf("MassA", _daughterMasses.first, "RelBreitWigner" + suffix);
  tr->createLeaf("MassB", _daughterMasses.second, "RelBreitWigner" + suffix);
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
  size_t check_nDouble = 15;
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
                       "Number of DoubleParameters does not match: " +
                       std::to_string(nDouble) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (nComplex != check_nComplex)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(nComplex) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (nMInteger != check_nMDouble)
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
  unsigned int spin = (unsigned int)paras.doubleParameter(2)->value();
  double d = paras.doubleParameter(3)->value();
  formFactorType ffType = formFactorType(paras.doubleParameter(4)->value());
  double ma = paras.doubleParameter(5)->value();
  double mb = paras.doubleParameter(6)->value();

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

void RelativisticBreitWigner::GetParameters(ParameterList &list) {
  AbstractDynamicalFunction::GetParameters(list);

  // We check of for each parameter if a parameter of the same name exists in
  // list. If so we check if both are equal and set the local parameter to the
  // parameter from the list. In this way we connect parameters that occur on
  // different positions in the amplitude.
  std::shared_ptr<DoubleParameter> tmp, width, radius;
  width = GetWidthParameter();
  radius = GetMesonRadiusParameter();
  try { // catch BadParameter
    tmp = FindParameter(width->name(), list);
    // catch and throw std::runtime_error due to failed parameter comparisson
    try {
      if (*tmp == *width)
        SetWidthParameter(tmp);
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.addParameter(width);
  }

  try { // catch BadParameter
    tmp = FindParameter(radius->name(), list);
    // catch and throw std::runtime_error due to failed parameter comparisson
    try {
      if (*tmp == *radius)
        SetMesonRadiusParameter(tmp);
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.addParameter(radius);
  }
}

void RelativisticBreitWigner::updateParameters(const ParameterList &list) {

  // Try to update mesonRadius
  std::shared_ptr<DoubleParameter> rad;
  try {
    rad = FindParameter(_mesonRadius->name(), list);
  } catch (std::exception &ex) {
  }
  if (rad)
    _mesonRadius->updateParameter(rad);

  // Try to update width
  std::shared_ptr<DoubleParameter> width;
  try {
    width = FindParameter(_width->name(), list);
  } catch (std::exception &ex) {
  }
  if (width)
    _width->updateParameter(width);

  return;
}
