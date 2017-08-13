// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
#include <numeric>
#include <iterator>
#include <vector>

#include "boost/property_tree/ptree.hpp"

#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/DecayDynamics/RelativisticBreitWigner.hpp"
#include "Physics/DecayDynamics/Coupling.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

std::shared_ptr<AbstractDynamicalFunction>
RelativisticBreitWigner::Factory(std::shared_ptr<PartList> partL,
                                 const boost::property_tree::ptree &pt) {
  auto obj = std::make_shared<RelativisticBreitWigner>();

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  LOG(trace) << "RelativisticBreitWigner::Factory() | Construction of " << name
             << ".";
  obj->SetName(name);
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
      auto width = ComPWA::DoubleParameterFactory(v.second);
      obj->SetWidthParameter(std::make_shared<DoubleParameter>(width));
    } else if (type == "MesonRadius") {
      auto mesonRadius = ComPWA::DoubleParameterFactory(v.second);
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
      << partProp.GetName() << " -> " << daughterNames.first << " + "
      << daughterNames.second;

  return std::static_pointer_cast<AbstractDynamicalFunction>(obj);
}

std::complex<double> RelativisticBreitWigner::Evaluate(const dataPoint &point,
                                                       int pos) const {
  std::complex<double> result = dynamicalFunction(
      point.GetValue(pos), _mass->GetValue(), _daughterMasses.first,
      _daughterMasses.second, _width->GetValue(), (double)_spin,
      _mesonRadius->GetValue(), _ffType);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return result;
}

bool RelativisticBreitWigner::CheckModified() const {
  if (AbstractDynamicalFunction::CheckModified())
    return true;
  if (_width->GetValue() != _current_width ||
      _mesonRadius->GetValue() != _current_mesonRadius) {
    SetModified();
    const_cast<double &>(_current_width) = _width->GetValue();
    const_cast<double &>(_current_mesonRadius) = _mesonRadius->GetValue();
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

std::shared_ptr<FunctionTree>
RelativisticBreitWigner::GetTree(const ParameterList &sample, int pos,
                                 std::string suffix) {

  //  int sampleSize = sample.GetMultiDouble(0)->GetNValues();

  std::shared_ptr<FunctionTree> tr(new FunctionTree());

  tr->createHead("RelBreitWigner" + suffix,
                 std::shared_ptr<Strategy>(new BreitWignerStrategy("")));

  tr->createLeaf("Mass", _mass, "RelBreitWigner" + suffix);
  tr->createLeaf("Width", _width, "RelBreitWigner" + suffix);
  tr->createLeaf("Spin", (double)_spin, "RelBreitWigner" + suffix);
  tr->createLeaf("MesonRadius", _mesonRadius, "RelBreitWigner" + suffix);
  tr->createLeaf("FormFactorType", _ffType, "RelBreitWigner" + suffix);
  tr->createLeaf("MassA", _daughterMasses.first, "RelBreitWigner" + suffix);
  tr->createLeaf("MassB", _daughterMasses.second, "RelBreitWigner" + suffix);
  tr->createLeaf("Data_mSq[" + std::to_string(pos) + "]",
                 sample.GetMultiDouble(pos), "RelBreitWigner" + suffix);

  return tr;
};

bool BreitWignerStrategy::execute(ParameterList &paras,
                                  std::shared_ptr<AbsParameter> &out) {
#ifndef NDEBUG
  // Check parameter type
  if (checkType != out->type())
    throw(WrongParType("BreitWignerStrat::execute() | "
                       "Output parameter is of type " +
                       std::string(ParNames[out->type()]) +
                       " and conflicts with expected type " +
                       std::string(ParNames[checkType])));

  // How many parameters do we expect?
  int check_nBool = 0;
  int check_nInt = 0;
  int check_nComplex = 0;
  int check_nDouble = 7;
  int check_nMDouble = 1;
  int check_nMComplex = 0;

  // Check size of parameter list
  if (paras.GetNBool() != check_nBool)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of BoolParameters does not match: " +
                       std::to_string(paras.GetNBool()) + " given but " +
                       std::to_string(check_nBool) + " expected."));
  if (paras.GetNInteger() != check_nInt)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of IntParameters does not match: " +
                       std::to_string(paras.GetNInteger()) + " given but " +
                       std::to_string(check_nInt) + " expected."));
  if (paras.GetNDouble() != check_nDouble)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of DoubleParameters does not match: " +
                       std::to_string(paras.GetNDouble()) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (paras.GetNComplex() != check_nComplex)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(paras.GetNComplex()) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (paras.GetNMultiDouble() != check_nMDouble)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of MultiDoubles does not match: " +
                       std::to_string(paras.GetNMultiDouble()) + " given but " +
                       std::to_string(check_nMDouble) + " expected."));
  if (paras.GetNMultiComplex() != check_nMComplex)
    throw(BadParameter("BreitWignerStrat::execute() | "
                       "Number of MultiComplexes does not match: " +
                       std::to_string(paras.GetNMultiComplex()) +
                       " given but " + std::to_string(check_nMComplex) +
                       " expected."));
#endif

  /** Get parameters from ParameterList:
   * We use the same order of the parameters as was used during tree
   * construction
   */
  double m0 = paras.GetDoubleParameter(0)->GetValue();
  double Gamma0 = paras.GetDoubleParameter(1)->GetValue();
  unsigned int spin = (unsigned int)paras.GetDoubleParameter(2)->GetValue();
  double d = paras.GetDoubleParameter(3)->GetValue();
  formFactorType ffType =
      formFactorType(paras.GetDoubleParameter(4)->GetValue());
  double ma = paras.GetDoubleParameter(5)->GetValue();
  double mb = paras.GetDoubleParameter(6)->GetValue();

  std::vector<double> mp = paras.GetMultiDouble(0)->GetValues();

  std::vector<std::complex<double>> results(mp.size(),
                                            std::complex<double>(0., 0.));

  // calc function for each point
  for (unsigned int ele = 0; ele < mp.size(); ele++) {
    try {
      results.at(ele) = RelativisticBreitWigner::dynamicalFunction(
          mp.at(ele), m0, ma, mb, Gamma0, spin, d, ffType);
    } catch (std::exception &ex) {
      LOG(error) << "BreitWignerStrategy::execute() | " << ex.what();
      throw(std::runtime_error("BreitWignerStrategy::execute() | "
                               "Evaluation of dynamic function failed!"));
    }
  }
  out =
      std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(), results));
  return true;
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
    tmp = list.GetDoubleParameter(width->GetName());
    // catch and throw std::runtime_error due to failed parameter comparisson
    try { 
      if (*tmp == *width)
        SetWidthParameter(tmp);
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.AddParameter(width);
  }

  try { // catch BadParameter
    tmp = list.GetDoubleParameter(radius->GetName());
    // catch and throw std::runtime_error due to failed parameter comparisson
    try {
      if (*tmp == *radius)
        SetMesonRadiusParameter(tmp);
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.AddParameter(radius);
  }
}

} /* namespace DecayDynamics */
} /* namespace Physics */
} /* namespace ComPWA */
