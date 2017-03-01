//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding flatte type resonance, removing root
//dependence
//-------------------------------------------------------------------------------
#include <numeric>
#include <algorithm>

#include "Core/PhysConst.hpp"
#include "Core/Functions.hpp"

// Boost header files go here
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/version.hpp>

#include "Physics/AmplitudeSum/NonResonant.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

#include "gsl/gsl_monte_vegas.h"

#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

using namespace boost::property_tree;

AmpSumIntensity::AmpSumIntensity(std::string name, normStyle ns,
                                 std::shared_ptr<Efficiency> eff,
                                 unsigned int nCalls)
    : AmpIntensity(name, eff), _maxFcnVal(0.), _calcMaxFcnVal(0), _normStyle(ns),
      _nCalls(nCalls) {
  return;
}

//! Copy constructor
AmpSumIntensity::AmpSumIntensity(const AmpSumIntensity &copy)
    : _maxFcnVal(copy._maxFcnVal), _calcMaxFcnVal(copy._calcMaxFcnVal),
      _normStyle(copy._normStyle), _nCalls(copy._nCalls) {
  // Copy efficiency
  _eff = copy._eff;

  // Deep copy of resonances
//  auto it = copy._ampList.begin();
//  for (; it != copy._ampList.end(); ++it) {
//    _ampList.push_back(std::shared_ptr<Resonance>((*it)->Clone()));
//  }

  // copy parameter list, but ensure that parameters are not added twice
//  int size = copy.params.GetNDouble();
//  for (unsigned int i = 0; i < size; ++i) {
//    bool found = 0;
//    std::string name = copy.params.GetDoubleParameter(i)->GetName();
//    for (unsigned int j = 0; j < params.GetNDouble(); ++j) {
//      if (params.GetDoubleParameter(j)->GetName() == name)
//        found = 1;
//    }
//    if (!found)
//      params.AddParameter(copy.params.GetDoubleParameter(i));
//  }
//  std::shared_ptr<DoubleParameter> r = copy.result.GetDoubleParameter(0);
//  result.AddParameter(
//      std::shared_ptr<DoubleParameter>(new DoubleParameter(*r)));

  // Check if memory addresses are different
  //	std::cout<<(result.GetDoubleParameter(0).get())<<"
  //"<<r.get()<<std::endl;
}

//! Clone function
AmpSumIntensity *AmpSumIntensity::Clone(std::string newName) const {
  auto tmp = (new AmpSumIntensity(*this));
  if (newName != "")
    tmp->SetName(newName);
  return tmp;
}

//==========================================================
//===================== OPERATORS ==========================
//==========================================================
/** Operator for coherent addition of amplitudes
 *
 * @param other
 * @return
 */
const AmpSumIntensity AmpSumIntensity::
operator+(const AmpSumIntensity &other) const {
  AmpSumIntensity ret(*this); // Make a copy of myself.
  ret += other;               // Use += to add other to the copy.
  return ret;
}

/** Operator for coherent addition of amplitudes
 *
 * @param rhs
 * @return
 */
AmpSumIntensity &AmpSumIntensity::operator+=(const AmpSumIntensity &rhs) {
  _name = _name + " + " + rhs._name;
  _ampList.insert(_ampList.end(), rhs._ampList.begin(), rhs._ampList.end());
  //    	params.insert(params.end(),
  //    rhs.params.begin(),rhs.params.begin());
  _calcMaxFcnVal = 0;
  if (_nCalls < rhs._nCalls)
    _nCalls = rhs._nCalls;
  return *this;
}

//==========================================================
//==================== PRINTING ============================
//==========================================================
void AmpSumIntensity::to_str() {
  std::stringstream outStr;
  outStr << "AmpSumIntensity: Printing resonances:\n";
//  auto it = GetResonanceItrFirst();
//  for (; it != GetResonanceItrLast(); ++it)
//    outStr << (*it)->to_str();

  LOG(info) << outStr.str();
  return;
}

void AmpSumIntensity::printFractions() {
  std::stringstream outStr;
  outStr << "Fit fractions for all amplitudes: \n";
  double sumFrac = 0;
//  auto it = GetResonanceItrFirst();
//  double norm = 1 / integral(GetResonanceItrList(), 0, _nCalls);
//  for (; it != GetResonanceItrLast(); ++it) {
//    double frac = (*it)->GetMagnitude() * norm;
//    sumFrac += frac;
//    outStr << std::setw(10) << (*it)->GetName() << ":    " << frac << "\n";
//  }

  outStr << std::setw(10) << " "
         << "    ==========\n";
  outStr << std::setw(10) << " "
         << "     " << sumFrac;
  LOG(info) << outStr.str();
  return;
}

double AmpSumIntensity::averageWidth() {
  double avWidth = 0;
  double sum = 0;
  for (int i = 0; i < _ampList.size(); i++) {
//    avWidth +=
//        std::norm(_ampList.at(i)->GetMagnitude()) * _ampList.at(i)->GetWidth();
//    sum += std::norm(_ampList.at(i)->GetMagnitude());
  }
  avWidth /= sum;
  return avWidth;
}

void AmpSumIntensity::SetPrefactor(std::complex<double> pre) {
  LOG(info) << "AmpSumIntensity::SetPreFactor() | "
               "Setting prefactor "
            << pre << " for all resonance of amplitude " << GetName() << "!";

  auto it = _ampList.begin();
  for (; it != _ampList.end(); ++it)
    (*it)->SetPreFactor(pre);
}

double AmpSumIntensity::GetMaximum(std::shared_ptr<Generator> gen) {
  if (!_calcMaxFcnVal)
    calcMaxVal(gen);
  return _maxFcnVal;
}

void AmpSumIntensity::calcMaxVal(std::shared_ptr<Generator> gen) {
  auto kin = dynamic_cast<DalitzKinematics *>(Kinematics::instance());

  double maxM23 = -999;
  double maxM13 = -999;
  double maxVal = 0;
  for (unsigned int i = 0; i < _nCalls; i++) {
    auto m13sq_limit = kin->GetMinMax(1);
    auto m23sq_limit = kin->GetMinMax(0);

    double m23sq =
        gen->getUniform() * (m23sq_limit.second - m23sq_limit.first) +
        m23sq_limit.first;
    double m13sq =
        gen->getUniform() * (m13sq_limit.second - m13sq_limit.first) +
        m13sq_limit.first;
    dataPoint point;
    try {
      Kinematics::instance()->FillDataPoint(1, 0, m13sq, m23sq, point);
    } catch (BeyondPhsp &ex) {
      if (i > 0)
        i--;
      continue;
    }
    double intens = Intensity(point);
    if (intens > maxVal) {
      maxM23 = m23sq;
      maxM13 = m13sq;
      maxVal = intens;
    }
  }
  _maxFcnVal = maxVal;
  _calcMaxFcnVal = 1;
  LOG(info) << "AmpSumIntensity::calcMaxVal() | "
               "calculated maximum of amplitude: "
            << _maxFcnVal << " at m23sq=" << maxM23 << "/m13sq=" << maxM13;
  return;
}

//=============================================================
//==================== NORMALIZATION/INTEGRATION ==============
//=============================================================
double AmpSumIntensity::Integral() const {
//  return Integral(GetResonanceItrList());
}

double AmpSumIntensity::Integral(std::vector<ampItr> _ampList) const{
  return AmpSumIntensity::integral(_ampList,
                                   0, // efficiency not included
                                   _nCalls);
}

double AmpSumIntensity::GetNormalization() const {
  double res=1;
//  double res = AmpSumIntensity::integral(GetResonanceItrList(), GetEfficiency(),
//                                         _nCalls);

  // check for NaN
  if (res != res)
    throw std::runtime_error("AmpSumIntensity::normalization() |"
                             "Result of amplitude " +
                             GetName() + " is NaN!");
  // check for inf
  if (std::isinf(res))
    throw std::runtime_error("AmpSumIntensity::normalization() |"
                             "Result of amplitude " +
                             GetName() + " is inf!");
  // check for zero
  if (res == 0)
    throw std::runtime_error("AmpSumIntensity::normalization() |"
                             "Result of amplitude " +
                             GetName() + " is zero!");

  return res;
}

struct GSLOpt_integral {
  std::vector<ampItr> resList;
  std::shared_ptr<Efficiency> eff;
};

double GSLWrapper_integral(double *x, size_t dim, void *param) {
  /* Calculation amplitude normalization (including efficiency).
   * We can not used (invmass, angle) for integration since the full point
   * has to be calculated.*/

  if (dim != 2)
    return 0;

  dataPoint point;

  try {
    Kinematics::instance()->FillDataPoint(0, 1, x[0], x[1], point);
  } catch (BeyondPhsp &ex) {
    return 0;
  }

  GSLOpt_integral *opt = static_cast<GSLOpt_integral *>(param);
  std::vector<ampItr> resList = opt->resList;
  std::shared_ptr<Efficiency> eff = opt->eff;

  std::complex<double> sum(0, 0);
  auto it = resList.begin();
  for (; it != resList.end(); ++it) {
    std::complex<double> val = (*(*it))->Evaluate(point);
    sum += val;
  }
  double pointEff = 1.0;
  if (eff)
    pointEff = eff->Evaluate(point);
  return std::norm(sum) * pointEff;
}

double AmpSumIntensity::integral(std::vector<ampItr> resList,
                                 std::shared_ptr<Efficiency> eff, int nCalls) {
  GSLOpt_integral par;
  par.resList = resList;
  par.eff = eff;

  /* Integration functionality was tested with a model with only one
   * normalized amplitude. The integration result is equal to the
   * amplitude coefficient^2.
   */
  size_t dim = 2;
  double res = 0.0, err = 0.0;

  auto kin = dynamic_cast<DalitzKinematics*>( Kinematics::instance() );

  // Set limits
  auto var1_limit = kin->GetMinMax(0);
  auto var2_limit = kin->GetMinMax(1);
  double xLimit_low[2] = {var1_limit.first, var2_limit.first};
  double xLimit_high[2] = {var1_limit.second, var2_limit.second};

  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default; // type of random generator
  gsl_rng *r = gsl_rng_alloc(T);           // random generator

  gsl_monte_function G = {&GSLWrapper_integral, dim, &par};

  /*	Choosing vegas algorithm here, because it is the most accurate:
   * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
   * 		 this should be sufficiency for most applications
   */
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
  gsl_monte_vegas_integrate(&G, xLimit_low, xLimit_high, dim, nCalls, r, s,
                            &res, &err);
  gsl_monte_vegas_free(s);

  // check for NaN
  if (res != res)
    throw std::runtime_error("AmpSumIntensity::integral() |"
                             "Result is NaN!");
  // check for inf
  if (std::isinf(res))
    throw std::runtime_error("AmpSumIntensity::integral() |"
                             "Result is inf!");

  LOG(debug) << "AmpSumIntensity::integrate() | Integration result"
                " for amplitude sum: "
             << res << "+-" << err << " relAcc [%]: " << 100 * err / res;

  return res;
}

double GSLWrapper_intfIntegral(double *x, size_t dim, void *param) {
  /* Calculation amplitude integral (excluding efficiency) */
  if (dim != 2)
    return 0;

  dataPoint point;
  try {
    Kinematics::instance()->FillDataPoint(0, 1, x[1], x[0], point);
  } catch (BeyondPhsp &ex) {
    return 0;
  }
  GSLOpt_integral *opt = static_cast<GSLOpt_integral *>(param);
  std::vector<ampItr> resList = opt->resList;
  std::shared_ptr<Efficiency> eff = opt->eff;

  double rm = 0;
  std::complex<double> sum(0, 0);
  auto it = resList.begin();
  for (; it != resList.end(); ++it) {
    std::complex<double> val = (*(*it))->Evaluate(point);
    rm += std::norm(val);
    sum += val;
  }
  double pointEff = 1.0;
  if (eff)
    pointEff = eff->Evaluate(point);
  return (std::norm(sum) - rm) * pointEff;
}

double
AmpSumIntensity::GetIntegralInterference(std::vector<ampItr> resList,
                                         unsigned int nCalls) {
  GSLOpt_integral par;
  par.resList = resList;

  /* Integration functionality was tested with a model with only one
   * normalized amplitude. The integration result is equal to the
   * amplitude coefficient^2.
   */
  size_t dim = 2;
  double res = 0.0, err = 0.0;

  // set limits: we assume that x[0]=m13sq and x[1]=m23sq
  auto kin = dynamic_cast<DalitzKinematics *>( Kinematics::instance() );

  // Set limits
  auto var1_limit = kin->GetMinMax(0);
  auto var2_limit = kin->GetMinMax(1);
  double xLimit_low[2] = {var1_limit.first, var2_limit.first};
  double xLimit_high[2] = {var1_limit.second, var2_limit.second};

  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default; // type of random generator
  gsl_rng *r = gsl_rng_alloc(T);           // random generator
  gsl_monte_function G = {&GSLWrapper_intfIntegral, dim, &par};

  /*	Choosing vegas algorithm here, because it is the most accurate:
   * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
   * 		 this should be sufficiency for most applications
   */
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
  gsl_monte_vegas_integrate(&G, xLimit_low, xLimit_high, 2, nCalls, r, s, &res,
                            &err);
  gsl_monte_vegas_free(s);

  return res;
}

double AmpSumIntensity::GetIntegralInterference(ampItr A,
                                                      ampItr B) const {
  std::vector<ampItr> par;
  par.push_back(A);
  par.push_back(B);
  double val = GetIntegralInterference(par, _nCalls);
  if (A == B)
    val /= 2;
  return val;
}

//=======================================================
//==================== EVALUATION =======================
//=======================================================
std::complex<double> AmpSumIntensity::Evaluate(const dataPoint &point) const {
  std::complex<double> AMPpdf(0, 0);
  auto it = _ampList.begin();
  for (; it != _ampList.end(); ++it) {
    try {
      AMPpdf += (*it)->Evaluate(point);
    } catch (std::exception &ex) {
      LOG(error) << "AmpSumIntensity::intensityNoEff() | "
                    "Failed to evaluate resonance "
                 << (*it)->GetName() << ": " << ex.what();
      throw;
    }
  }
  return AMPpdf;
}

double AmpSumIntensity::IntensityNoEff(const dataPoint &point) const{
  return std::norm(Evaluate(point));
}

double AmpSumIntensity::Intensity(const dataPoint &point) const {
  double ampNoEff;
  try {
    ampNoEff = IntensityNoEff(point);
  } catch (BadParameter &ex) {
    LOG(error) << "AmpSumIntensity::Intensity() | Can not "
                  "obtain parameter from ParameterList 'result'!";
    throw;
  }
  return ampNoEff * _eff->Evaluate(point);
}

double
AmpSumIntensity::sliceIntensity(dataPoint &dataP, ParameterList &par,
                                std::complex<double> *reso, unsigned int nResos,
                                double N, unsigned int nF0, unsigned int nF2) {
  double AMPpdf = 0;
  // TODO: implement slice fit
  //	if(Kinematics::instance()->isWithinPhsp(dataP))
  //	AMPpdf = totAmp.evaluateSlice(dataP, reso, nResos,5, N, nF0, nF2);
  if (AMPpdf != AMPpdf) {
    LOG(error) << "Error AmpSumIntensity: Intensity is not a number!!";
    AMPpdf = 0;
  }
  double eff = _eff->Evaluate(dataP);
  return AMPpdf * eff;
}

//============== FIT FRACTIONS ================
ampItr findResonancePartner(std::shared_ptr<AmpSumIntensity> amp,
                                  ampItr res) {
  auto name = (*res)->GetName();
  auto it = amp->begin();
  for (; it != amp->end(); ++it) { // fill matrix
    if (it == res)
      continue;
    auto name2 = (*it)->GetName();
    if (name2.find(name) != std::string::npos)
      return it;
  }
  return res;
}

void AmpSumIntensity::GetFitFractions(ParameterList &parList) {
  GetFitFractions(parList, this);
}

void AmpSumIntensity::GetFitFractions(ParameterList &parList,
                                      const AmpSumIntensity *ampConst) {
  // Work around: Remove const
  AmpSumIntensity* amp = const_cast<AmpSumIntensity*>(ampConst);

  std::string ampName = amp->GetName();

  // Select only resonace without "_CP" in name for calculation of fit fraction
  auto it = amp->begin();
  std::vector<ampItr> _ampList;
  for (; it != amp->begin(); ++it) {
    if ((*it)->GetName().find("_CP") != std::string::npos)
      continue;
    _ampList.push_back(it);
  }
  if (!_ampList.size())
    throw std::runtime_error(
        "AmpSumIntensity::GetFitFractions() | "
        "No resonance are selected for calculation if fit fractions!");

  /* Unbinned efficiency correction in the FunctionTree does not provide
   * an integral w/o efficiency correction. We have to calculate it here. */
  double norm = 1.0;
  try {
    norm = amp->Integral(_ampList);
    // norm = amp->GetIntegral();
  } catch (std::exception &ex) {
    LOG(error) << "FitResult::calcFraction() | "
                  "Normalization can't be calculated: "
               << ex.what();
    throw;
  }

  LOG(debug) << "FitResult::calcFraction() | "
                "Amplitude "
             << ampName << " Norm=" << norm;

  // Start loop over resonances
  auto itit = _ampList.begin();
  for (; itit != _ampList.end(); ++itit) {
    // Current resonance iterator
    ampItr cReso = (*itit);

    // We search for a partner resonance and add it to the integral
    //		auto it2 = findResonancePartner(amp, it);

    // GetIntegralInterference returns the integal Int( A*B+B*A ),
    // including the complex coefficienct
    //		double nom = amp->GetIntegralInterference(it,it);
    //		if( it != it2 ){// Int |A+B|^2 = |A|^2 + |B|^2 + A*B + B*A
    //			double tmp22 = amp->GetIntegralInterference(it2,it2);
    //			double tmp12 = amp->GetIntegralInterference(it,it2);
    //			LOG(debug) << "FitResult::calcFraction() | Calculating"
    //					<<" amplitude integral for composed amplitudes
    //"
    //					<<(*it)->GetName()<<" and "<<(*it2)->GetName()<<":
    //"
    //					<<"(11) "<<nom <<" (22) "<<tmp22 <<" (12)
    //"<<tmp12
    //					<<" Total: "<<nom+tmp22+tmp12;
    //			nom += tmp22;
    //			nom += tmp12;
    //		} else {
    //		LOG(debug) << "FitResult::calcFraction() | Resonance "
    //				"integal for "<<(*it)->GetName()<<": "<<nom;
    //		}
    std::vector<ampItr> thisAmp;
    thisAmp.push_back(cReso);

    // Calculate resonance integral. This includes the magnitude^2.
//    double nom = amp->Integral(thisAmp);
    double nom = 1.0;

    LOG(debug) << "FitResult::calcFraction() | Resonance "
                  "integal for "
               << (*cReso)->GetName() << ": " << nom;

    std::string resName = ampName + " " + (*cReso)->GetName() + "_FF";
    std::shared_ptr<DoubleParameter> magPar = (*cReso)->GetMagnitude();
    double mag = magPar->GetValue(); // value of magnitude
    double magError = 0;
    if (magPar->HasError())
      magError = magPar->GetError(); // error of magnitude

    parList.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter(
        resName, nom / norm, std::fabs(2 * (nom / mag) / norm * magError))));
  }
}

//=========================================================
//================== ACCESS to resonances =================
//=========================================================

int AmpSumIntensity::GetIdOfResonance(std::string name) {
  for (unsigned int i = 0; i < _ampList.size(); i++)
    if (_ampList.at(i)->GetName() == name)
      return i;
  return -999;
}

std::string AmpSumIntensity::GetNameOfResonance(unsigned int id) {
  if (id > _ampList.size())
    throw std::runtime_error("AmpSumIntensity::GetNameOfResonance() | "
                             "Invalid resonance ID=" +
                             std::to_string(id) + "! Resonance not found?");
  return _ampList.at(id)->GetName();
}

std::shared_ptr<Amplitude> AmpSumIntensity::GetAmplitude(std::string name) {
  int id = GetIdOfResonance(name);
  return GetAmplitude(id);
}

std::shared_ptr<Amplitude> AmpSumIntensity::GetAmplitude(unsigned int id) {
  if (id > _ampList.size())
    throw std::runtime_error("AmpSumIntensity::GetResonance() | "
                             "Invalid resonance ID=" +
                             std::to_string(id) + "! Resonance not found?");
  return _ampList.at(id);
}

//=====================================
//========== FunctionTree =============
//=====================================
//! Getter function for function tree
std::shared_ptr<FunctionTree>
AmpSumIntensity::GetTree(ParameterList &sample, ParameterList &phspSample,
                         ParameterList &toySample) {
  unsigned int effId = Kinematics::instance()->GetNVars();
  unsigned int weightId = Kinematics::instance()->GetNVars() + 1;

  int phspSampleSize = phspSample.GetMultiDouble(0)->GetNValues();
  /* We assume that the total efficiency of the phsp variable is stored in
   * each event. This does not influence the result but a correct efficiency
   * given consistent results with the GSL integration. */
  //	double phspSampleEff = phspSample.GetMultiDouble(effId)->GetValue(0);

  std::shared_ptr<MultiDouble> weightPhsp = phspSample.GetMultiDouble(weightId);
  double sumWeights =
      std::accumulate(weightPhsp->Begin(), weightPhsp->End(), 0.0);
  std::shared_ptr<MultiDouble> eff = phspSample.GetMultiDouble(effId);

  std::shared_ptr<Strategy> mmultStrat(new MultAll(ParType::MCOMPLEX));
  std::shared_ptr<Strategy> mmultDStrat(new MultAll(ParType::MDOUBLE));
  std::shared_ptr<Strategy> multiDoubleAddStrat(new AddAll(ParType::MDOUBLE));
  std::shared_ptr<Strategy> multiComplexAddStrat(new AddAll(ParType::MCOMPLEX));
  std::shared_ptr<Strategy> msqStrat(new AbsSquare(ParType::MDOUBLE));
  std::shared_ptr<Strategy> mlogStrat(new LogOf(ParType::MDOUBLE));
  std::shared_ptr<Strategy> multStrat(new MultAll(ParType::COMPLEX));
  std::shared_ptr<Strategy> multDStrat(new MultAll(ParType::DOUBLE));
  std::shared_ptr<Strategy> addStrat(new AddAll(ParType::DOUBLE));
  std::shared_ptr<Strategy> addComplexStrat(new AddAll(ParType::COMPLEX));
  std::shared_ptr<Strategy> sqStrat(new AbsSquare(ParType::DOUBLE));
  std::shared_ptr<Strategy> logStrat(new LogOf(ParType::DOUBLE));
  std::shared_ptr<Strategy> complStrat(new Complexify(ParType::COMPLEX));
  std::shared_ptr<Strategy> invStrat(new Inverse(ParType::DOUBLE));

  //------------Setup Tree---------------------
  std::shared_ptr<FunctionTree> tr(new FunctionTree());

  tr->createHead(GetName(), mmultDStrat);
  tr->createNode("AmpSq", msqStrat, GetName());
  tr->insertTree(setupBasicTree(sample, toySample), "AmpSq");

  // Normalization
  tr->createNode("N", invStrat, GetName()); // 1/normLH
  // normLH = phspVolume/N_{mc} |T_{evPHSP}|^2
  tr->createNode("normFactor", multDStrat, "N");
  // sumAmp = \sum_{evPHSP} |T_{evPHSP}|^2
  tr->createNode("sumAmp", addStrat, "normFactor");
  tr->createLeaf("phspVolume", Kinematics::instance()->GetPhspVolume(),
                 "normFactor");
  tr->createLeaf("InvNmc", 1 / ((double)sumWeights), "normFactor");
  tr->createNode("IntensPhspEff", mmultDStrat, "sumAmp", phspSampleSize,
                 false);                       //|T_{ev}|^2
  tr->createLeaf("eff", eff, "IntensPhspEff"); // efficiency
  tr->createLeaf("weightPhsp", weightPhsp, "IntensPhspEff");
  tr->createNode("IntensPhsp", msqStrat, "IntensPhspEff", phspSampleSize,
                 false); //|T_{ev}|^2
  tr->insertTree(setupBasicTree(phspSample, toySample, "_norm"), "IntensPhsp");

  return tr;
}

std::shared_ptr<FunctionTree>
AmpSumIntensity::setupBasicTree(ParameterList &sample,
                                ParameterList &phspSample, std::string suffix) {
  int sampleSize = sample.GetMultiDouble(0)->GetNValues();
  int phspSampleSize = phspSample.GetMultiDouble(0)->GetNValues();

  if (sampleSize == 0) {
    LOG(error) << "AmpSumIntensity::setupBasicTree() | "
                  "Data sample empty!";
    return std::shared_ptr<FunctionTree>();
  }
  if (phspSampleSize == 0) {
    LOG(error) << "AmpSumIntensity::setupBasicTree() | "
                  "Phsp sample empty!";
    return std::shared_ptr<FunctionTree>();
  }

  //------------Setup Tree---------------------
  std::shared_ptr<FunctionTree> newTree(new FunctionTree());

  //----Strategies needed
  std::shared_ptr<AddAll> maddStrat(new AddAll(ParType::MCOMPLEX));

  newTree->createHead("Amplitude" + suffix, maddStrat, sampleSize);

  auto it = _ampList.begin();
  for (; it != _ampList.end(); ++it) {
    std::shared_ptr<FunctionTree> resTree =
        (*it)->GetTree(sample, phspSample, phspSample,"" + (*it)->GetName());
    if (!resTree->sanityCheck())
      throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
                               "Resonance tree didn't pass sanity check!");
    resTree->recalculate();
    newTree->insertTree(resTree, "Amplitude" + suffix);
  }

  LOG(debug) << "AmpSumIntensity::setupBasicTree(): tree constructed!!";
  return newTree;
}
} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */
