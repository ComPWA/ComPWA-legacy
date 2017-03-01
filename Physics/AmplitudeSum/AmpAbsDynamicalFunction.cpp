//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------

#include <stdlib.h>
#include "gsl/gsl_monte.h"
#include "gsl/gsl_monte_plain.h"
#include "gsl/gsl_monte_miser.h"
#include "gsl/gsl_monte_vegas.h"

#include "Core/PhysConst.hpp"
#include "Physics/AmplitudeSum/DalitzKinematics.hpp"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(normStyle nS, int calls)
    : Amplitude(""), _modified(1), _nCalls(calls), _normStyle(nS),
      _ffType(formFactorType::BlattWeisskopf), _parity(+1), _cparity(0) {}

AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(
    const char *name, unsigned int varIdA, unsigned int varIdB,
    std::shared_ptr<DoubleParameter> mag,
    std::shared_ptr<DoubleParameter> phase,
    std::shared_ptr<DoubleParameter> mass, Spin spin, Spin m, Spin n, int P,
    int C, std::string mother, std::string particleA, std::string particleB,
    std::shared_ptr<DoubleParameter> mesonR,  //  meson radius
    std::shared_ptr<DoubleParameter> motherR, //  mother radius
    formFactorType type, int nCalls, normStyle nS)
    : Amplitude(name), _modified(1), _nCalls(nCalls),
      _nameMother(mother), _name1(particleA), _name2(particleB), _mag(mag),
      _phase(phase), _mass(mass), _normStyle(nS), _ffType(type),
      _mesonRadius(mesonR), _motherRadius(motherR), _subSys(varIdA),
      _spin(spin), _m(m), _n(n), _parity(P), _cparity(C),
      _wignerD(varIdB, spin) {
  initialize();
}

AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(
    const char *name, unsigned int varIdA, unsigned int varIdB,
    std::shared_ptr<DoubleParameter> mag,
    std::shared_ptr<DoubleParameter> phase,
    std::shared_ptr<DoubleParameter> mass, Spin spin, Spin m, Spin n, int P,
    int C, std::string mother, std::string particleA, std::string particleB,
    formFactorType type, int nCalls, normStyle nS)
    : Amplitude(name), _modified(1), _nCalls(nCalls),
      _nameMother(mother), _name1(particleA), _name2(particleB), _mag(mag),
      _phase(phase), _mass(mass), _normStyle(nS), _ffType(type),
      _mesonRadius(std::make_shared<DoubleParameter>(name, 1.0)),
      _motherRadius(std::make_shared<DoubleParameter>(name, 1.0)),
      _subSys(varIdA), _spin(spin), _m(m), _n(n), _parity(P), _cparity(C),
      _wignerD(varIdB, spin) {
  initialize();
}

std::string AmpAbsDynamicalFunction::to_str() const {
  std::stringstream str;
  str << "AmpAbsDynamicalFunction | " << _name
      << " nCalls=" << _nCalls << " varId1=" << GetVarIdA()
      << " varId2=" << GetVarIdB() << std::endl
      << " J=" << (double)_spin << " P=" << _parity << " C=" << _cparity
      << " ffType=" << _ffType << " mother: " << _nameMother
      << " particleA: " << _name1 << " particleB: " << _name2 << std::endl;
  str << " normStyle=" << _normStyle << " modified?" << _modified
      << " Prefactor=" << _preFactor << std::endl;
  str << "Parameters:" << std::endl;
  str << _mag->to_str() << std::endl;
  str << _phase->to_str() << std::endl;
  str << _mass->to_str() << std::endl;
  str << _mesonRadius->to_str() << std::endl;
  str << _motherRadius->to_str() << std::endl;

  return str.str();
}

void AmpAbsDynamicalFunction::initialize() {
  try {
    _M = PhysConst::Instance()->FindParticle(_nameMother).GetMass();
  } catch (...) {
    throw BadConfig("AmpAbsDynamicalFunction::initialize() | "
                    "Can not obtain mass of mother particle: " +
                    _nameMother);
  }

  try {
    _mass1 = PhysConst::Instance()->FindParticle(_name1).GetMass();
  } catch (...) {
    throw BadConfig("AmpAbsDynamicalFunction::initialize() | "
                    "Can not obtain mass of daughter 1: " +
                    _name1);
  }
  try {
    _mass2 = PhysConst::Instance()->FindParticle(_name2).GetMass();
  } catch (...) {
    throw BadConfig("AmpAbsDynamicalFunction::initialize() | "
                    "DCan not obtain mass of daughter 2: " +
                    _name2);
  }
}

AmpAbsDynamicalFunction::~AmpAbsDynamicalFunction() {}


  void AmpAbsDynamicalFunction::CheckModified() const
  {
if(_mass->GetValue() != tmp_mass){
  const_cast<bool&>(_modified) = 1;
		const_cast<double&>(tmp_mass) = _mass->GetValue();
}
  }
  
std::complex<double> AmpAbsDynamicalFunction::Evaluate(const dataPoint &point) const {
  CheckModified();
  std::complex<double> res(0, 0);
  try {
    res = EvaluateAmp(point);
  } catch (std::exception &ex) {
    LOG(error) << "AmpAbsDynamicalFunction::Evaluate() | "
                  "Failed to evaluate dynamic part of resonance "
               << GetName() << " at"
                               " point:"
               << std::endl
               << point << std::endl
               << ex.what();
    throw;
  }
  double ang = 0;
  try {
    ang = EvaluateWignerD(point);
  } catch (std::exception &ex) {
    LOG(error) << "AmpAbsDynamicalFunction::Evaluate() | "
                  "Failed to evaluate angular part of resonance "
               << GetName() << " at"
                               " point:"
               << std::endl
               << point << std::endl
               << ex.what();
    throw;
  }
  res = (GetPreFactor() * std::polar( GetMagnitudeValue(), GetPhaseValue() ) * GetNormalization() * res * ang);

#ifndef NDEBUG
  // check for NaN
  if (std::isnan(res.real()) || std::isnan(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::Evaluate() | Result of"
                             " resonance " +
                             GetName() + " is NaN!");
  // check for inf
  if (std::isinf(res.real()) || std::isinf(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::Evaluate() | Result of"
                             " resonance " +
                             GetName() + " is inf!");
#endif

  return res;
}

double evalAmp(double *x, size_t dim, void *param) {
  /* We need a wrapper here because a eval() is a member function of
   * AmpAbsDynamicalFunction and can therefore not be referenced. But
   * gsl_monte_function expects a function reference. As third parameter
   * we pass the reference to the current instance of AmpAbsDynamicalFunction
   */
  if (dim != 2)
    return 0;

  auto amp = static_cast<AmpAbsDynamicalFunction *>(param);
  dataPoint point;

  try {
    Kinematics::instance()->FillDataPoint(0, 1, x[0], x[1], point);
  } catch (BeyondPhsp &ex) {
    return 0;
  }

  //		int idA = amp->GetVarIdA();
  //		int idB = amp->GetVarIdB();
  //		if( !Kinematics::instance()->IsWithinBoxPhsp(idA, idB, x[0], x[1])
  //)
  //			return 0;
  //		point.setVal(idA, x[0]);
  //		point.setVal(idB, x[1]);

  std::complex<double> res(0, 0);
  try {
    res = amp->EvaluateAmp(point);
  } catch (std::exception &ex) {
    LOG(error) << "AmpAbsDynamicalFunction -> evalAmp() | "
                  "Amplitude can not be evaluated at point "
               << point << "! " << ex.what();
    throw;
  }
  // include angular distribution in normalization
  res *= amp->EvaluateWignerD(point);

  return (std::norm(res)); // integrate over |F|^2
}

double AmpAbsDynamicalFunction::Integral() const {
  size_t dim = 2;
  double res = 0.0, err = 0.0;

  auto kin = dynamic_cast<DalitzKinematics*>( Kinematics::instance() );

  //	auto var1_limit = kin->GetMinMax( GetVarIdA() );
  //	auto var2_limit = kin->GetMinMax( GetVarIdB() );
  //	double vol = (var1_limit.second-var1_limit.first)
  //			*(var2_limit.second-var2_limit.first);
  auto var1_limit = kin->GetMinMax(0);
  auto var2_limit = kin->GetMinMax(1);
  //	double vol = kin->GetPhspVolume();
  double vol = 1.0;
  double xLimit_low[2] = {var1_limit.first, var2_limit.first};
  double xLimit_high[2] = {var1_limit.second, var2_limit.second};

  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default; // type of random generator
  gsl_rng *r = gsl_rng_alloc(T);           // random generator
  gsl_monte_function F = {&evalAmp, dim,
                          const_cast<AmpAbsDynamicalFunction *>(this)};

  // Test function: result should be 1
  // gsl_monte_function F = {&twoDimGaussian,dim, new int()};

  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
  gsl_monte_vegas_integrate(&F, xLimit_low, xLimit_high, 2, _nCalls, r, s, &res,
                            &err);
  gsl_monte_vegas_free(s);

#ifndef NDEBUG
  // check for NaN
  if (std::isnan(res))
    throw std::runtime_error("AmpAbsDynamicalFunction::integral() |"
                             "Result of resonance " +
                             GetName() + " is NaN!");
  // check for inf
  if (std::isinf(res))
    throw std::runtime_error("AmpAbsDynamicalFunction::integral() |"
                             "Result of resonance " +
                             GetName() + " is inf!");
  // check for zero
  if (res == 0)
    throw std::runtime_error("AmpAbsDynamicalFunction::integral() |"
                             "Result of resonance " +
                             GetName() + " is zero!");

  LOG(debug) << "AmpAbsDynamicalFunction::integral() | "
                "Integration result for |"
             << _name << "|^2: " << res << "+-" << err
             << " relAcc [%]: " << 100 * err / res;
#endif

  return res / vol;
}

double AmpAbsDynamicalFunction::GetNormalization() const {
  if (_normStyle == normStyle::none)
    return 1.0;
  double norm = 1 / sqrt(Integral());

#ifndef NDEBUG
  // check for NaN
  if (std::isnan(norm))
    throw std::runtime_error("AmpAbsDynamicalFunction::GetNormalization() |"
                             "Result of resonance " +
                             GetName() + " is NaN!");
  // check for inf
  if (std::isinf(norm))
    throw std::runtime_error("AmpAbsDynamicalFunction::GetNormalization() |"
                             "Result of resonance " +
                             GetName() + " is inf!");
  // check for zero
  if (norm == 0)
    throw std::runtime_error("AmpAbsDynamicalFunction::GetNormalization() |"
                             "Result of resonance " +
                             GetName() + " is zero!");
#endif

  return norm;
}

double AmpAbsDynamicalFunction::GetTotalIntegral() const {
  // TODO: add test case to assure that the integral is one
  if (_normStyle == normStyle::one)
    return std::norm(GetPreFactor());
  return TotalIntegral();
}

double eval(double *x, size_t dim, void *param) {
  /* We need a wrapper here because evaluate() is a member function of
   * AmpAbsDynamicalFunction
   * and can therefore not be referenced. But gsl_monte_function expects a
   * function reference.
   * As third parameter we pass the reference to the current instance of
   * AmpAbsDynamicalFunction
   */
  if (dim != 2)
    return 0;

  auto amp = static_cast<AmpAbsDynamicalFunction *>(param);
  dataPoint point;

  try {
    Kinematics::instance()->FillDataPoint(0, 1, x[0], x[1], point);
  } catch (BeyondPhsp &ex) {
    return 0;
  }

  //	int idA = amp->GetVarIdA();
  //	int idB = amp->GetVarIdB();
  //	if( !Kinematics::instance()->IsWithinBoxPhsp(idA, idB, x[0], x[1]) )
  //		return 0;
  //	point.setVal(idA, x[0]);
  //	point.setVal(idB, x[1]);

  //	std::complex<double> res = amp->EvaluateAmp(point);
  //	double ang = amp->EvaluateWignerD(point);
  //	double norm = amp->GetNormalization();
  //	return ( std::norm(res*ang*norm) ); //integrate over |F|^2
  return (std::norm(amp->Evaluate(point) /
                    std::polar( amp->GetMagnitudeValue(), amp->GetPhaseValue() )) ); // integrate over |F|^2
}

double AmpAbsDynamicalFunction::TotalIntegral() const {
  size_t dim = 2;
  double res = 0.0, err = 0.0;

   auto kin = dynamic_cast<DalitzKinematics *>( Kinematics::instance() );

  auto var1_limit = kin->GetMinMax(0);
  auto var2_limit = kin->GetMinMax(1);
  //	auto var1_limit = kin->GetMinMax( GetVarIdA() );
  //	auto var2_limit = kin->GetMinMax( GetVarIdB() );
  double xLimit_low[2] = {var1_limit.first, var2_limit.first};
  double xLimit_high[2] = {var1_limit.second, var2_limit.second};

  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default; // type of random generator
  gsl_rng *r = gsl_rng_alloc(T);           // random generator
  gsl_monte_function F = {&eval, dim,
                          const_cast<AmpAbsDynamicalFunction *>(this)};

  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
  gsl_monte_vegas_integrate(&F, xLimit_low, xLimit_high, 2, _nCalls, r, s, &res,
                            &err);
  gsl_monte_vegas_free(s);

#ifndef NDEBUG
  // check for NaN
  if (std::isnan(res))
    throw std::runtime_error("AmpAbsDynamicalFunction::totalIntegral() |"
                             "Result of resonance " +
                             GetName() + " is NaN!");
  // check for inf
  if (std::isinf(res))
    throw std::runtime_error("AmpAbsDynamicalFunction::totalIntegral() |"
                             "Result of resonance " +
                             GetName() + " is inf!");
  // check for zero
  if (res == 0)
    throw std::runtime_error("AmpAbsDynamicalFunction::totalIntegral() |"
                             "Result of resonance " +
                             GetName() + " is zero!");
#endif

  return res;
}

std::complex<double> AmpAbsDynamicalFunction::widthToCoupling(
    double mSq, double mR, double width, double ma, double mb, double spin,
    double mesonRadius, formFactorType type) {
  double sqrtS = sqrt(mSq);

  // calculate gammaA(s_R)
  std::complex<double> gammaA(1, 0); // spin==0
  if (spin > 0) {
    std::complex<double> qValue = Kinematics::qValue(mR, ma, mb);
    double ffR =
        Kinematics::FormFactor(mR, ma, mb, spin, mesonRadius, qValue, type);
    std::complex<double> qR = std::pow(qValue, spin);
    gammaA = ffR * qR;
  }

  // calculate phsp factor
  std::complex<double> rho = Kinematics::phspFactor(sqrtS, ma, mb);

  std::complex<double> denom = gammaA * sqrt(rho);
  std::complex<double> res = std::complex<double>(sqrt(mR * width), 0) / denom;

#ifndef NDEBUG
  // check for NaN
  //	if( std::isnan(res.real()) || std::isnan(res.imag()) )
  //		throw std::runtime_error("AmpAbsDynamicalFunction::widthToCoupling()
  //| "
  //				"Result is NaN!");
  // check for inf
  if (std::isinf(res.real()) || std::isinf(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::widthToCoupling() | "
                             "Result is inf!");
#endif

  return res;
}

std::complex<double> AmpAbsDynamicalFunction::couplingToWidth(
    double mSq, double mR, double g, double ma, double mb, double spin,
    double mesonRadius, formFactorType type) {
  double sqrtM = sqrt(mSq);
  std::complex<double> phspFactor = Kinematics::phspFactor(sqrtM, ma, mb);

  return couplingToWidth(mSq, mR, g, ma, mb, spin, mesonRadius, type,
                         phspFactor);
}

std::complex<double> AmpAbsDynamicalFunction::couplingToWidth(
    double mSq, double mR, double g, double ma, double mb, double spin,
    double mesonRadius, formFactorType type, std::complex<double> phspFactor) {
  // calculate gammaA(s_R)
  std::complex<double> gammaA(1, 0); // spin==0
  if (spin > 0 || type == formFactorType::CrystalBarrel) {
    std::complex<double> qValue = Kinematics::qValue(mR, ma, mb);
    double ffR =
        Kinematics::FormFactor(mR, ma, mb, spin, mesonRadius, qValue, type);
    std::complex<double> qR = std::pow(qValue, spin);
    gammaA = ffR * qR;
  }
  // calculate phsp factor
  std::complex<double> res = std::norm(gammaA) * g * g * phspFactor / mR;

#ifndef NDEBUG
  // check for NaN
  if (std::isnan(res.real()) || std::isnan(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::couplingToWidth() | "
                             "Result is NaN!");
  // check for inf
  if (std::isinf(res.real()) || std::isinf(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::couplingToWidth() | "
                             "Result is inf!");
#endif

  return res;
}

//_____________________________________________________________________________

double twoDimGaussian(double *z, size_t dim, void *param) {
  if (dim != 2)
    return 0;
  /* test environment for numeric integration:
   * 	Calculating integral of normalized gaussian:
   * 	f(x,y) = A exp( - (x-x0)^2/(2*sigmaX^2) + (y-y0)^2/(2*sigmaY^2)
   * 	with A=1/(2*pi*sigmaX*sigmaY) this function is normalized to 1
   */
  double x = z[0];
  double y = z[1];
  // mean and width need to be adjusted according to final state kinematics
  double x0 = 1.1, y0 = 1.1;           // mean
  double sigmaX = 0.01, sigmaY = 0.01; // width
  double pi = ComPWA::PhysConst::Instance()->FindConstant("Pi").GetValue();

  double result = exp(-(x - x0) * (x - x0) / (2 * sigmaX * sigmaX) -
                      (y - y0) * (y - y0) / (2 * sigmaY * sigmaY));
  result /= 2 * pi * sigmaY * sigmaX;
  return result;
}

//_____________________________________________________________________________
std::shared_ptr<FunctionTree> couplingToWidthStrat::SetupTree(
    std::shared_ptr<MultiDouble> mSq, std::shared_ptr<DoubleParameter> mR,
    std::shared_ptr<DoubleParameter> g, double ma, double mb, Spin spin,
    std::shared_ptr<DoubleParameter> mesonRadius, formFactorType type,
    std::string suffix) {

  std::shared_ptr<couplingToWidthStrat> thisStrat(new couplingToWidthStrat);

  std::string stratName = "couplingToWidth" + suffix;
  //------------Setup Tree---------------------
  std::shared_ptr<FunctionTree> newTree(new FunctionTree());

  newTree->createHead(stratName, thisStrat, mSq->GetNValues());
  newTree->createLeaf("mass", mR, stratName);
  newTree->createLeaf("g", g, stratName);
  newTree->createLeaf("massA", ma, stratName);
  newTree->createLeaf("massB", mb, stratName);
  newTree->createLeaf("spin", (double)spin, stratName);
  newTree->createLeaf("mesonRadius", mesonRadius, stratName);
  newTree->createLeaf("ffType", (double)type, stratName);

  newTree->insertTree(phspFactorStrat::SetupTree(mSq, ma, mb), stratName);
  newTree->createLeaf("mSq", mSq, stratName);

  return newTree;
}

bool couplingToWidthStrat::execute(ParameterList &paras,
                                   std::shared_ptr<AbsParameter> &out) {
#ifndef NDEBUG
  // Check parameter type
  if (checkType != out->type())
    throw(WrongParType("couplingToWidthStrat::execute() | "
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
  int check_nMComplex = 1;

  // Check size of parameter list
  if (paras.GetNBool() != check_nBool)
    throw(BadParameter("couplingToWidthStrat::execute() | "
                       "Number of BoolParameters does not match: " +
                       std::to_string(paras.GetNBool()) + " given but " +
                       std::to_string(check_nBool) + " expected."));
  if (paras.GetNInteger() != check_nInt)
    throw(BadParameter("couplingToWidthStrat::execute() | "
                       "Number of IntParameters does not match: " +
                       std::to_string(paras.GetNInteger()) + " given but " +
                       std::to_string(check_nInt) + " expected."));
  if (paras.GetNDouble() != check_nDouble)
    throw(BadParameter("couplingToWidthStrat::execute() | "
                       "Number of DoubleParameters does not match: " +
                       std::to_string(paras.GetNDouble()) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (paras.GetNComplex() != check_nComplex)
    throw(BadParameter("couplingToWidthStrat::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(paras.GetNComplex()) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (paras.GetNMultiDouble() != check_nMDouble)
    throw(BadParameter("couplingToWidthStrat::execute() | "
                       "Number of MultiDoubles does not match: " +
                       std::to_string(paras.GetNMultiDouble()) + " given but " +
                       std::to_string(check_nMDouble) + " expected."));
  if (paras.GetNMultiComplex() != check_nMComplex)
    throw(BadParameter("couplingToWidthStrat::execute() | "
                       "Number of MultiComplexes does not match: " +
                       std::to_string(paras.GetNMultiComplex()) +
                       " given but " + std::to_string(check_nMComplex) +
                       " expected."));
#endif

  /* Get parameters from ParameterList:
   * We use the same order of the parameters as was used during tree
   * construction */
  double mR = paras.GetDoubleParameter(0)->GetValue();
  double g = paras.GetDoubleParameter(1)->GetValue();
  double ma = paras.GetDoubleParameter(2)->GetValue();
  double mb = paras.GetDoubleParameter(3)->GetValue();
  unsigned int spin = paras.GetDoubleParameter(4)->GetValue();
  double mesonRadius = paras.GetDoubleParameter(5)->GetValue();
  formFactorType ffType =
      formFactorType(paras.GetDoubleParameter(6)->GetValue());

  std::vector<double> mp = paras.GetMultiDouble(0)->GetValues();
  std::vector<std::complex<double>> phspFactors =
      paras.GetMultiComplex(0)->GetValues();

  std::vector<std::complex<double>> results(mp.size(),
                                            std::complex<double>(0., 0.));
  // calc function for each point
  for (unsigned int ele = 0; ele < mp.size(); ele++) {
    try {
      results.at(ele) = AmpAbsDynamicalFunction::couplingToWidth(
          mp.at(ele), mR, g, ma, mb, spin, mesonRadius, ffType,
          phspFactors.at(ele));
    } catch (std::exception &ex) {
      LOG(error) << "couplingToWidthStrat::execute() | " << ex.what();
      throw(std::runtime_error("couplingToWidthStrat::execute() | "
                               "Evaluation of dynamic function failed!"));
    }
  }
  out =
      std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(), results));
  return true;
}
} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */
