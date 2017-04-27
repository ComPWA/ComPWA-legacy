//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//     Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include <numeric>

#include "boost/property_tree/ptree.hpp"

#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/HelicityFormalism/RelativisticBreitWigner.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::complex<double> RelativisticBreitWigner::Evaluate(const dataPoint &point,
                                                       int pos) const {
  return EvaluateNoNorm(point.GetValue(pos)) * GetNormalization();
}

std::complex<double> RelativisticBreitWigner::EvaluateNoNorm(double mSq) const {
  std::complex<double> result = dynamicalFunction(
      mSq, _mass->GetValue(), _massA, _massB, _width->GetValue(), (double)_spin,
      _mesonRadius->GetValue(), _ffType);
  return result;
}

double RelativisticBreitWigner::GetNormalization() const {
  CheckModified();
  if (GetModified()) {
    const_cast<double &>(_current_integral) =
        AbstractDynamicalFunction::Integral();
    SetModified(false);
  }
  return 1. / _current_integral;
}

void RelativisticBreitWigner::CheckModified() const {
  AbstractDynamicalFunction::CheckModified();
  if (_width->GetValue() != _current_width ||
      _mesonRadius->GetValue() != _current_mesonRadius) {
    SetModified();
    const_cast<double &>(_current_width) = _width->GetValue();
    const_cast<double &>(_current_mesonRadius) = _mesonRadius->GetValue();
  }
  return;
}

std::complex<double> RelativisticBreitWigner::dynamicalFunction(
    double mSq, double mR, double ma, double mb, double width, unsigned int J,
    double mesonRadius, formFactorType ffType) {

  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  auto phspFactorSqrtS = Kinematics::phspFactor(sqrtS, ma, mb);
  auto phspFactormR = Kinematics::phspFactor(mR, ma, mb);
  
  // Check if we have an event which is exactly at the phase space boundary
  if(phspFactorSqrtS == std::complex<double>(0,0))
    return std::complex<double>(0,0);
  
  std::complex<double> qTerm = std::pow(( phspFactorSqrtS / phspFactormR ) * mR / sqrtS,
                                        (2 * J + 1));
  double barrier =
      HelicityKinematics::FormFactor(sqrtS, ma, mb, J, mesonRadius, ffType) /
      HelicityKinematics::FormFactor(mR, ma, mb, J, mesonRadius, ffType);


  // Calculate coupling constant to final state
  std::complex<double> g_final =
      widthToCoupling(mSq, mR, width, ma, mb, J, mesonRadius, ffType);

  /*Coupling constant from production reaction. In case of a particle decay
   * the production coupling doesn't depend in energy since the CM energy
   * is in the (RC) system fixed to the mass of the decaying particle */
  double g_production = 1;

  std::complex<double> denom = std::complex<double>(mR * mR - mSq, 0) +
                               (-1.0) * i * sqrtS * (width * qTerm * barrier);

  std::complex<double> result = g_final * g_production / denom;

  assert(!std::isnan(result.real()));
  assert(!std::isnan(result.imag()));

  return result;
}

std::shared_ptr<RelativisticBreitWigner>
RelativisticBreitWigner::Factory(const boost::property_tree::ptree &pt) {
  LOG(trace) << "RelativisticBreitWigner::Factory() | Construction....";
  auto obj = std::make_shared<RelativisticBreitWigner>();

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  obj->SetName(name);
  auto partProp = PhysConst::Instance()->FindParticle(name);
  obj->SetMass(std::make_shared<DoubleParameter>(partProp.GetMassPar()));

  auto decayTr = partProp.GetDecayInfo();
  if (partProp.GetDecayType() != "relativisticBreitWigner")
    throw std::runtime_error(
        "RelativisticBreitWigner::Factory() | Decay type does not match! ");

  auto spin = partProp.GetSpin();
  obj->SetSpin(spin);

  auto ffType = formFactorType(decayTr.get<int>("FormFactor.<xmlattr>.type"));
  obj->SetFormFactorType(ffType);

  auto width = ComPWA::DoubleParameterFactory(decayTr.get_child("Width"));
  obj->SetWidth(std::make_shared<DoubleParameter>(width));
  auto mesonRadius =
      ComPWA::DoubleParameterFactory(decayTr.get_child("MesonRadius"));
  obj->SetMesonRadius(std::make_shared<DoubleParameter>(mesonRadius));

  // Get masses of decay products
  auto decayProducts = pt.get_child("DecayProducts");
  std::vector<std::string> daughterNames;
  for (auto i : decayProducts) {
    daughterNames.push_back(i.second.get<std::string>("<xmlattr>.Name"));
  }
  if (daughterNames.size() != 2)
    throw boost::property_tree::ptree_error(
        "AmpWignerD::Factory() | Expect exactly two decay products (" +
        std::to_string(decayProducts.size()) + " given)!");

  obj->SetDecayMassA(
      PhysConst::Instance()->FindParticle(daughterNames.at(0)).GetMass());
  obj->SetDecayMassB(
      PhysConst::Instance()->FindParticle(daughterNames.at(1)).GetMass());
  obj->SetDecayNameA(daughterNames.at(0));
  obj->SetDecayNameB(daughterNames.at(1));

  LOG(trace)
      << "RelativisticBreitWigner::Factory() | Construction of the decay "
      << partProp.GetName() << " -> " << daughterNames.at(0) << " + "
      << daughterNames.at(1);

  return obj;
}

/**! Setup function tree */
std::shared_ptr<FunctionTree>
RelativisticBreitWigner::GetTree(ParameterList &sample,
                                 ParameterList &toySample, int pos,
                                 std::string suffix) {

  int sampleSize = sample.GetMultiDouble(0)->GetNValues();
  int phspSampleSize = toySample.GetMultiDouble(0)->GetNValues();

  std::shared_ptr<FunctionTree> tr(new FunctionTree());

  tr->createHead("DynamicalFunction",
                 std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)));

  tr->createNode("RelBreitWigner",
                 std::shared_ptr<Strategy>(new BreitWignerStrategy("")),
                 "DynamicalFunction", sampleSize);
  tr->createLeaf("Mass", _mass, "RelBreitWigner");               // m0
  tr->createLeaf("Width", _width, "RelBreitWigner");             // resWidth
  tr->createLeaf("Spin", (double)_spin, "RelBreitWigner");       // spin
  tr->createLeaf("MesonRadius", _mesonRadius, "RelBreitWigner"); // d
  tr->createLeaf("FormFactorType", _ffType, "RelBreitWigner");   // d
  tr->createLeaf("MassA", _massA, "RelBreitWigner");             // ma
  tr->createLeaf("MassB", _massB, "RelBreitWigner");             // mb
  tr->createLeaf("Data_mSq[" + std::to_string(pos) + "]",
                 sample.GetMultiDouble(pos), "RelBreitWigner");

  // Normalization
  tr->createNode("Normalization",
                 std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)),
                 "DynamicalFunction"); // 1/normLH
  tr->createNode("Integral",
                 std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
                 "Normalization");
  tr->createNode("Sum", std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                 "Integral");
  tr->createLeaf("Range", _limits.second - _limits.first, "Integral");
  tr->createLeaf("InvNmc", 1 / ((double)phspSampleSize), "Integral");
  tr->createNode("IntensPhsp",
                 std::shared_ptr<Strategy>(new AbsSquare(ParType::MDOUBLE)),
                 "Sum", phspSampleSize,
                 false); //|T_{ev}|^2
  tr->createNode("NormRelBreitWigner",
                 std::shared_ptr<Strategy>(new BreitWignerStrategy("")),
                 "IntensPhsp", phspSampleSize);
  tr->createLeaf("Mass", _mass, "NormRelBreitWigner");               
  tr->createLeaf("Width", _width, "NormRelBreitWigner");
  tr->createLeaf("Spin", (double)_spin, "NormRelBreitWigner");
  tr->createLeaf("MesonRadius", _mesonRadius, "NormRelBreitWigner"); 
  tr->createLeaf("FormFactorType", _ffType, "NormRelBreitWigner");   
  tr->createLeaf("MassA", _massA, "NormRelBreitWigner");             
  tr->createLeaf("MassB", _massB, "NormRelBreitWigner");
  tr->createLeaf("PhspSample_mSq[" + std::to_string(pos) + "]",
                 toySample.GetMultiDouble(pos),
                 "NormRelBreitWigner"); // mc

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

  /* We check of for each parameter if a parameter of the same name exists in
   * list. If so we check if both are equal and set the local parameter to the
   * parameter from the list. In this way we connect parameters that occur on
   * different positions in the amplitude.
   */
  std::shared_ptr<DoubleParameter> tmp, width, radius;
  width = GetWidth();
  radius = GetMesonRadius();
  try { // catch BadParameter
    tmp = list.GetDoubleParameter(width->GetName());
    try { // catch and throw std::runtime_error due to failed parameter
          // comparisson
      if (*tmp == *width)
        SetWidth(tmp);
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.AddParameter(width);
  }

  try { // catch BadParameter
    tmp = list.GetDoubleParameter(radius->GetName());
    try { // catch and throw std::runtime_error due to failed parameter
          // comparisson
      if (*tmp == *radius)
        SetMesonRadius(tmp);
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.AddParameter(radius);
  }
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
