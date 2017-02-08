 
    
//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding correct g1s
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include <math.h>
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

using Physics::DPKinematics::DalitzKinematics;

AmpFlatteRes::AmpFlatteRes(normStyle nS, int calls)
    : AmpAbsDynamicalFunction(nS, calls) {}

AmpFlatteRes::AmpFlatteRes(
    const char *name, unsigned int varIdA, unsigned int varIdB,
    std::shared_ptr<DoubleParameter> mag,
    std::shared_ptr<DoubleParameter> phase,
    std::shared_ptr<DoubleParameter> mass, Spin spin, Spin m, Spin n, int P,
    int C, std::string mother, std::string particleA, std::string particleB,
    std::shared_ptr<DoubleParameter> mesonRadius,
    std::shared_ptr<DoubleParameter> motherRadius,
    std::shared_ptr<DoubleParameter> g1, std::shared_ptr<DoubleParameter> g2,
    std::string g2_idA, std::string g2_idB, std::shared_ptr<DoubleParameter> g3,
    std::string g3_idA, std::string g3_idB, formFactorType type, int nCalls,
    normStyle nS)
    : AmpAbsDynamicalFunction(name, varIdA, varIdB, mag, phase, mass, spin, m,
                              n, P, C, mother, particleA, particleB,
                              mesonRadius, motherRadius, type, nCalls, nS),
      _g1(g1), _g2(g2), _g3(g3), _g2_idA(g2_idA), _g2_idB(g2_idB),
      _g3_idA(g3_idA), _g3_idB(g3_idB)

{
  tmp_g1 = _g1->GetValue();
  tmp_g2 = _g2->GetValue();
  tmp_g3 = _g3->GetValue();

  initialize();
}

AmpFlatteRes::~AmpFlatteRes() {}

void AmpFlatteRes::initialize() {
  AmpAbsDynamicalFunction::initialize();

  try {
    _g2_massA = PhysConst::instance()->findParticle(_g2_idA).mass_;
    _g2_massB = PhysConst::instance()->findParticle(_g2_idB).mass_;
  } catch (std::exception &ex) {
    LOG(error) << "AmpFlatteRes::initialize() | "
                  "Masses of 2nd channel can not be set!";
    throw;
  }
  try {
    _g3_massA = PhysConst::instance()->findParticle(_g3_idA).mass_;
    _g3_massB = PhysConst::instance()->findParticle(_g3_idB).mass_;
  } catch (std::exception &ex) {
    // if g3 is not set we don't need to masses either
    if (_g3->GetValue()) {
      LOG(error) << "AmpFlatteRes::initialize() | "
                    "Masses of 3rd channel can not be set!";
      throw;
    }
  }
}

void AmpFlatteRes::Configure(boost::property_tree::ptree::value_type const &v,
                             ParameterList &list) {
  if (v.first != "Flatte")
    throw BadConfig("");

  boost::property_tree::ptree pt = v.second;
  AmpAbsDynamicalFunction::Configure(v, list);

  // Coupling to channel 1 (mandatory)
  auto tmp_g1_fix = pt.get<bool>("g1_fix", 1);
  auto tmp_g1_min = pt.get<double>("g1_min", 0.);
  auto tmp_g1_max = pt.get<double>("g1_max", 10.);
  auto tmp_g1_name = pt.get_optional<std::string>("g1_name");
  if (!tmp_g1_name) {
    auto tmp_g1 = pt.get_optional<double>("g1");
    if (!tmp_g1)
      throw BadParameter("AmpFlatteRes::Configure() | "
                         "g1 for " +
                         _name + " not specified!");
    _g1 = std::shared_ptr<DoubleParameter>(new DoubleParameter(
        "g1_" + _name, tmp_g1.get(), tmp_g1_min, tmp_g1_max));
    _g1->FixParameter(tmp_g1_fix);
    if (_enable)
      list.AddParameter(_g1);
    _g1_writeByName = 0;
  } else {
    try {
      _g1 = list.GetDoubleParameter(tmp_g1_name.get());
      _g1_writeByName = 1;
    } catch (BadParameter &ex) {
      if (!_enable) {
        _g1 = std::shared_ptr<DoubleParameter>(
            new DoubleParameter(tmp_g1_name.get(), 0.0));
        _g1_writeByName = 1;
      } else {
        LOG(error) << "AmpFlatteRes::Configure() | "
                      "Requesting parameter "
                   << tmp_g1_name.get() << " but"
                                           " was not found in parameter list. "
                                           "Quit since parameter is mandatory!";
        throw;
      }
    }
  }
  // Coupling to channel 2 (mandatory)
  auto tmp_g2_fix = pt.get<bool>("g2_fix", 1);
  auto tmp_g2_min = pt.get<double>("g2_min", 0.);
  auto tmp_g2_max = pt.get<double>("g2_max", 10.);
  auto tmp_g2_name = pt.get_optional<std::string>("g2_name");
  if (!tmp_g2_name) {
    auto tmp_g2 = pt.get_optional<double>("g2");
    if (!tmp_g2)
      throw BadParameter("AmpFlatteRes::Configure() | "
                         "g2 for " +
                         _name + " not specified!");
    _g2 = std::shared_ptr<DoubleParameter>(new DoubleParameter(
        "g2_" + _name, tmp_g2.get(), tmp_g2_min, tmp_g2_max));
    _g2->FixParameter(tmp_g2_fix);
    if (_enable)
      list.AddParameter(_g2);
    _g2_writeByName = 0;
  } else {
    try {
      _g2 = list.GetDoubleParameter(tmp_g2_name.get());
      _g2_writeByName = 1;
    } catch (BadParameter &ex) {
      if (!_enable) {
        _g2 = std::shared_ptr<DoubleParameter>(
            new DoubleParameter(tmp_g2_name.get(), 0.0));
        _g2_writeByName = 1;
      } else {
        LOG(error) << "AmpFlatteRes::Configure() | "
                      "Requesting parameter "
                   << tmp_g2_name.get() << " but"
                                           " was not found in parameter list. "
                                           "Quit since parameter is mandatory!";
        throw;
      }
    }
  }

  auto tmp_g2_massA = pt.get_optional<std::string>("g2_part1");
  if (!tmp_g2_massA)
    throw BadParameter("AmpFlatteRes::Configure() | "
                       "g2_massA for " +
                       _name + " not specified!");
  _g2_idA = tmp_g2_massA.get();

  auto tmp_g2_massB = pt.get_optional<std::string>("g2_part2");
  if (!tmp_g2_massB)
    throw BadParameter("AmpFlatteRes::Configure() | "
                       "g2_massB for " +
                       _name + " not specified!");
  _g2_idB = tmp_g2_massB.get();

  // Coupling to channel 3 (optional)
  auto tmp_g3_fix = pt.get<bool>("g3_fix", 1);
  auto tmp_g3_min = pt.get<double>("g3_min", 0.);
  auto tmp_g3_max = pt.get<double>("g3_max", 10.);
  auto tmp_g3_name = pt.get_optional<std::string>("g3_name");
  if (!tmp_g3_name) {
    auto tmp_g3 = pt.get<double>("g3", 0.0);
    _g3 = std::shared_ptr<DoubleParameter>(
        new DoubleParameter("g3_" + _name, tmp_g3, tmp_g3_min, tmp_g3_max));
    _g3->FixParameter(tmp_g3_fix);
    _g3_writeByName = 0;
    auto tmp_g3_massA = pt.get<std::string>("g3_part1", "");
    _g3_idA = tmp_g3_massA;
    auto tmp_g3_massB = pt.get<std::string>("g3_part2", "");
    _g3_idB = tmp_g3_massB;

    if (_enable && _g3->GetValue() != 0) {
      LOG(info) << "AmpFlatteRes::Configure() | Adding 3rd channel: " << _name
                << " parName=" << _g3->GetName() << " g3=" << _g3->GetValue()
                << " g3_massA=" << _g3_massA << " g3_massB=" << _g3_massB;
      list.AddParameter(_g3);
    }
  } else {
    try {
      _g3 = list.GetDoubleParameter(tmp_g3_name.get());
      _g3_writeByName = 1;
      auto tmp_g3_massA = pt.get<std::string>("g3_part1", "");
      _g3_idA = tmp_g3_massA;
      auto tmp_g3_massB = pt.get<std::string>("g3_part2", "");
      _g3_idB = tmp_g3_massB;
      LOG(info) << "AmpFlatteRes::Configure() | Adding 3rd channel:" << _name
                << " parName=" << _g3->GetName() << " g3=" << _g3->GetValue()
                << " g3_massA=" << _g3_massA << " g3_massB=" << _g3_massB;
    } catch (BadParameter &ex) {
      if (!_enable) {
        _g3 = std::shared_ptr<DoubleParameter>(
            new DoubleParameter(tmp_g3_name.get(), 0.0));
        _g3_writeByName = 1;
      } else {
        LOG(error) << "AmpFlatteRes::Configure() | "
                      "Requesting parameter "
                   << tmp_g3_name.get()
                   << " but"
                      " was not found in parameter list. "
                      "Continue since parameter is not mandatory!";
        _g3 = std::shared_ptr<DoubleParameter>(
            new DoubleParameter("g3_" + _name, 0.0));
      }
    }
  }

  initialize();

  return;
}

void AmpFlatteRes::Save(boost::property_tree::ptree &pt) {
  boost::property_tree::ptree amp;
  AmpAbsDynamicalFunction::put(amp);
  if (_g1_writeByName) {
    amp.put("g1_name", _g1->GetName());
  } else {
    amp.put("g1", _g1->GetValue());
    amp.put("g1_fix", _g1->IsFixed());
    amp.put("g1_min", _g1->GetMinValue());
    amp.put("g1_max", _g1->GetMaxValue());
  }
  if (_g2_writeByName) {
    amp.put("g2_name", _g2->GetName());
  } else {
    amp.put("g2", _g2->GetValue());
    amp.put("g2_fix", _g2->IsFixed());
    amp.put("g2_min", _g2->GetMinValue());
    amp.put("g2_max", _g2->GetMaxValue());
  }
  amp.put("g2_part1", _g2_idA);
  amp.put("g2_part2", _g2_idB);
  if (_g3->GetValue() != 0) {
    if (_g3_writeByName) {
      amp.put("g3_name", _g3->GetName());
    } else {
      amp.put("g3", _g3->GetValue());
      amp.put("g3_fix", _g3->IsFixed());
      amp.put("g3_min", _g3->GetMinValue());
      amp.put("g3_max", _g3->GetMaxValue());
    }
    amp.put("g3_part1", _g3_idA);
    amp.put("g3_part2", _g3_idB);
  }

  pt.add_child("Flatte", amp);

  return;
}

std::string AmpFlatteRes::to_str() const {
  std::string dynAmp = AmpAbsDynamicalFunction::to_str();
  std::stringstream str;
  str << _g1->to_str() << std::endl;
  str << _g2->to_str() << std::endl;
  if (_g3->GetValue())
    str << _g3->to_str() << std::endl;
  str << "massB1=" << _g2_massA << " massB2=" << _g2_massB;
  if (_g3->GetValue())
    str << " massC1=" << _g3_massA << " massC2=" << _g3_massB << std::endl;

  return dynAmp + str.str();
}

void AmpFlatteRes::CheckModified() {
  AmpAbsDynamicalFunction::CheckModified();
  if (_g1->GetValue() != tmp_g1 || _g2->GetValue() != tmp_g2 ||
      _g3->GetValue() != tmp_g3) {
    SetModified();
    tmp_g1 = _g1->GetValue();
    tmp_g2 = _g2->GetValue();
    tmp_g3 = _g3->GetValue();
  }
  return;
}

double AmpFlatteRes::GetIntegral() {
  CheckModified();
  if (_modified) {
    tmp_integral = integral();
    _modified = 0;
  }
  return tmp_integral;
}

std::complex<double> AmpFlatteRes::EvaluateAmp(dataPoint &point) {
  double mSq = point.getVal(_subSys);

  std::complex<double> result;
  try {
    result = dynamicalFunction(
        mSq, _mass->GetValue(), _mass1, _mass2, _g1->GetValue(), _g2_massA,
        _g2_massB, _g2->GetValue(), _g3_massA, _g3_massB, _g3->GetValue(),
        (double)_spin, _mesonRadius->GetValue(), _ffType);
  } catch (std::exception &ex) {
    LOG(error) << "AmpFlatteRes::EvaluateAmp() | "
                  "Dynamical function can not be evalutated: "
               << ex.what();
    throw;
  }
  return result;
}
std::complex<double> AmpFlatteRes::dynamicalFunction(
    double mSq, double mR, double gA, std::complex<double> termA,
    std::complex<double> termB, std::complex<double> termC) {
  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  /* Coupling constant from production reaction. In case of a particle decay
   * the production. coupling doesn't depend in energy since the CM energy
   * is in the (RC) system fixed to the mass of the decaying particle
   */
  double g_production = 1;

  std::complex<double> denom = std::complex<double>(mR * mR - mSq, 0) +
                               (-1.0) * i * sqrtS * (termA + termB + termC);

  std::complex<double> result =
      std::complex<double>(gA * g_production, 0) / denom;

#ifdef DEBUG
  if (std::isnan(result.real()) || std::isnan(result.imag())) {
    std::cout << "AmpFlatteRes::dynamicalFunction() | " << mR << " " << mSq
              << " " << termA << " " << termB << " " << termC << std::endl;
    return 0;
  }
#endif

  return result;
}

std::complex<double>
AmpFlatteRes::dynamicalFunction(double mSq, double mR, double massA1,
                                double massA2, double gA, double massB1,
                                double massB2, double couplingB, double massC1,
                                double massC2, double couplingC, unsigned int J,
                                double mesonRadius, formFactorType ffType) {
  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

// channel A - signal channel
std:
  complex<double> gammaA, qTermA, termA;
  double barrierA;
  // break-up momentum
  barrierA =
      Kinematics::FormFactor(sqrtS, massA1, massA2, J, mesonRadius, ffType) /
      Kinematics::FormFactor(mR, massA1, massA2, J, mesonRadius, ffType);
  // convert coupling to partial width of channel A
  gammaA = couplingToWidth(mSq, mR, gA, massA1, massA2, J, mesonRadius, ffType);
  // including the factor qTermA, as suggested by PDG, leads to an amplitude
  // that doesn't converge.
  //		qTermA = Kinematics::qValue(sqrtS,massA1,massA2) /
  //Kinematics::qValue(mR,massA1,massA2);
  qTermA = std::complex<double>(1, 0);
  termA = gammaA * barrierA * barrierA * std::pow(qTermA, (double)2 * J + 1);

  // channel B - hidden channel
  std::complex<double> gammaB, qTermB, termB;
  double barrierB, gB;
  // break-up momentum
  barrierB =
      Kinematics::FormFactor(sqrtS, massB1, massB2, J, mesonRadius, ffType) /
      Kinematics::FormFactor(mR, massB1, massB2, J, mesonRadius, ffType);
  gB = couplingB;
  // convert coupling to partial width of channel B
  gammaB = couplingToWidth(mSq, mR, gB, massB1, massB2, J, mesonRadius, ffType);
  //		qTermB = Kinematics::qValue(sqrtS,massB1,massB2) /
  //Kinematics::qValue(mR,massB1,massB2);
  qTermB = std::complex<double>(1, 0);
  termB = gammaB * barrierB * barrierB * std::pow(qTermB, (double)2 * J + 1);

  // channel C - hidden channel
  std::complex<double> gammaC, qTermC, termC;
  double barrierC, gC;
  if (couplingC != 0.0) {
    // break-up momentum
    barrierC =
        Kinematics::FormFactor(sqrtS, massC1, massC2, J, mesonRadius, ffType) /
        Kinematics::FormFactor(mR, massC1, massC2, J, mesonRadius, ffType);
    gC = couplingC;
    // convert coupling to partial width of channel C
    gammaC =
        couplingToWidth(mSq, mR, gC, massC1, massC2, J, mesonRadius, ffType);
    //		qTermC = Kinematics::qValue(sqrtS,massC1,massC2) /
    //Kinematics::qValue(mR,massC1,massC2);
    qTermC = std::complex<double>(1, 0);
    termC = gammaC * barrierC * barrierC * std::pow(qTermC, (double)2 * J + 1);
  }

  return dynamicalFunction(mSq, mR, gA, termA, termB, termC);
}

std::shared_ptr<FunctionTree> AmpFlatteRes::SetupTree(ParameterList &sample,
                                                      ParameterList &toySample,
                                                      std::string suffix) {
  DalitzKinematics *kin =
      dynamic_cast<DalitzKinematics *>(Kinematics::instance());
  double phspVol = kin->GetPhspVolume();

  int sampleSize = sample.GetMultiDouble(0)->GetNValues();
  int toySampleSize = toySample.GetMultiDouble(0)->GetNValues();

  LOG(info) << "AmpFlatteRes::setupBasicTree() | " << _name
            << " nEvents=" << sampleSize << " nPhspEvents=" << toySampleSize;

  //------------Setup Tree---------------------
  std::shared_ptr<FunctionTree> newTree(new FunctionTree());

  //----Strategies needed
  std::shared_ptr<MultAll> mmultStrat(new MultAll(ParType::MCOMPLEX));
  std::shared_ptr<AbsSquare> msqStrat(new AbsSquare(ParType::MDOUBLE));
  std::shared_ptr<MultAll> multDStrat(new MultAll(ParType::DOUBLE));
  std::shared_ptr<AddAll> addStrat(new AddAll(ParType::DOUBLE));
  std::shared_ptr<Complexify> complStrat(new Complexify(ParType::COMPLEX));
  std::shared_ptr<Inverse> invStrat(new Inverse(ParType::DOUBLE));
  std::shared_ptr<SquareRoot> sqRootStrat(new SquareRoot(ParType::DOUBLE));

  //----Add Nodes
  std::shared_ptr<FlatteStrategy> flatteStrat(new FlatteStrategy(_name));
  std::shared_ptr<WignerDStrategy> angdStrat(new WignerDStrategy(_name));

  // R = (mag,phase)*(Flatte)*(WignerD)*(Norm)
  newTree->createHead("Reso_" + _name, mmultStrat, sampleSize);
  newTree->createNode("PreFactor_" + _name, complStrat, "Reso_" + _name);
  newTree->createLeaf("IntensPre_" + _name, std::abs(_prefactor),
                      "PreFactor_" + _name);
  newTree->createLeaf("PhasePre_" + _name, std::arg(_prefactor),
                      "PreFactor_" + _name);
  newTree->createNode("C_" + _name, complStrat, "Reso_" + _name); // c
  newTree->createLeaf("Intens_" + _name, _mag, "C_" + _name);     // r
  newTree->createLeaf("Phase_" + _name, _phase, "C_" + _name);    // phi

  // Angular distribution
  if ((double)_spin)
    newTree->insertTree(_wignerD.SetupTree(sample, suffix), "Reso_" + _name);

  // Flatte
  newTree->insertTree(FlatteStrategy::SetupTree(
                          "FlatteRes_" + _name,
                          sample.GetMultiDouble(GetVarIdA()), _mass, _g1,
                          _mass1, _mass2, _g2, _g2_massA, _g2_massB, _g3,
                          _g3_massA, _g3_massB, _spin, _mesonRadius, _ffType),
                      "Reso_" + _name);

  // Normalization
  if (_normStyle == normStyle::none) {
    newTree->createLeaf("N_" + _name, 1., "Reso_" + _name);
  } else {
    newTree->createNode("N_" + _name, sqRootStrat, "Reso_" + _name);
    newTree->createNode("NSq_" + _name, multDStrat, "N_" + _name);
    newTree->createLeaf("PhspSize_" + _name, toySampleSize, "NSq_" + _name);
    newTree->createLeaf("PhspVolume_" + _name, 1 / phspVol, "NSq_" + _name);
    newTree->createNode("InvSum_" + _name, invStrat, "NSq_" + _name);
    newTree->createNode("Sum_" + _name, addStrat, "InvSum_" + _name);
    newTree->createNode("AbsVal_" + _name, msqStrat, "Sum_" + _name);

    newTree->createNode("NormReso_" + _name, mmultStrat, "AbsVal_" + _name,
                        toySampleSize); // BW

    // Angular distribution (Normalization)
    if ((double)_spin)
      newTree->insertTree(_wignerD.SetupTree(toySample, suffix),
                          "NormReso_" + _name);

    // Flatte (Normalization)
    newTree->insertTree(FlatteStrategy::SetupTree(
                            "NormFlatte_" + _name,
                            toySample.GetMultiDouble(GetVarIdA()), _mass, _g1,
                            _mass1, _mass2, _g2, _g2_massA, _g2_massB, _g3,
                            _g3_massA, _g3_massB, _spin, _mesonRadius, _ffType),
                        "NormReso_" + _name);
  }
  return newTree;
}

std::shared_ptr<FunctionTree> FlatteStrategy::SetupTree(
    std::string name, std::shared_ptr<MultiDouble> mSq,
    std::shared_ptr<DoubleParameter> mR, std::shared_ptr<DoubleParameter> g,
    double ma, double mb, std::shared_ptr<DoubleParameter> g2, double g2_ma,
    double g2_mb, std::shared_ptr<DoubleParameter> g3, double g3_ma,
    double g3_mb, Spin spin, std::shared_ptr<DoubleParameter> mesonRadius,
    formFactorType type) {
  std::shared_ptr<FlatteStrategy> thisStrat(new FlatteStrategy(name));

  std::string stratName = name;
  //------------Setup Tree---------------------
  std::shared_ptr<FunctionTree> newTree(new FunctionTree());

  newTree->createHead(stratName, thisStrat, mSq->GetNValues());

  newTree->createLeaf("mass", mR, stratName);
  newTree->createLeaf("g", g, stratName);

  // Partial width of channel A
  newTree->insertTree(couplingToWidthStrat::SetupTree(mSq, mR, g, ma, mb, spin,
                                                      mesonRadius, type, "_cA"),
                      stratName);
  newTree->insertTree(
      barrierStrat::SetupTree(mSq, mR, ma, mb, spin, mesonRadius, type, "_bA"),
      stratName);

  // Partial width of channel B
  newTree->insertTree(couplingToWidthStrat::SetupTree(mSq, mR, g2, g2_ma, g2_mb,
                                                      spin, mesonRadius, type,
                                                      "_cB"),
                      stratName);
  newTree->insertTree(barrierStrat::SetupTree(mSq, mR, g2_ma, g2_mb, spin,
                                              mesonRadius, type, "_bB"),
                      stratName);

  // Partial width of channel C (optional)
  if (g3->GetValue()) {
    newTree->insertTree(
        couplingToWidthStrat::SetupTree(mSq, mR, g3, g3_ma, g3_mb, spin,
                                        mesonRadius, type, "_cC"),
        stratName);
    newTree->insertTree(barrierStrat::SetupTree(mSq, mR, g3_ma, g3_mb, spin,
                                                mesonRadius, type, "_bC"),
                        stratName);
  } else {
    std::shared_ptr<MultiDouble> vecmd(
        new MultiDouble("zero", std::vector<double>(mSq->GetNValues(), 0.0)));
    std::shared_ptr<MultiComplex> veccd(new MultiComplex(
        "zero", std::vector<std::complex<double>>(mSq->GetNValues(),
                                                  std::complex<double>(0, 0))));
    newTree->createLeaf("gammaC", veccd, stratName);
    newTree->createLeaf("barrierSqC", vecmd, stratName);
  }

  // invariant masses
  newTree->createLeaf("mSq", mSq, stratName);

  return newTree;
}

bool FlatteStrategy::execute(ParameterList &paras,
                             std::shared_ptr<AbsParameter> &out) {
#ifdef DEBUG
  // Check parameter type
  if (checkType != out->type())
    throw(WrongParType("FlatteStrategy::execute() | "
                       "Output parameter is of type " +
                       std::string(ParNames[out->type()]) +
                       " and conflicts with expected type " +
                       std::string(ParNames[checkType])));

  // How many parameters do we expect?
  int check_nBool = 0;
  int check_nInt = 0;
  int check_nComplex = 0;
  int check_nDouble = 2;
  int check_nMDouble = 4;
  int check_nMComplex = 3;

  // Check size of parameter list
  if (paras.GetNBool() != check_nBool)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of BoolParameters does not match: " +
                       std::to_string(paras.GetNBool()) + " given but " +
                       std::to_string(check_nBool) + " expected."));
  if (paras.GetNInteger() != check_nInt)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of IntParameters does not match: " +
                       std::to_string(paras.GetNInteger()) + " given but " +
                       std::to_string(check_nInt) + " expected."));
  if (paras.GetNDouble() != check_nDouble)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of DoubleParameters does not match: " +
                       std::to_string(paras.GetNDouble()) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (paras.GetNComplex() != check_nComplex)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(paras.GetNComplex()) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (paras.GetNMultiDouble() != check_nMDouble)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of MultiDoubles does not match: " +
                       std::to_string(paras.GetNMultiDouble()) + " given but " +
                       std::to_string(check_nMDouble) + " expected."));
  if (paras.GetNMultiComplex() != check_nMComplex)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of MultiComplexes does not match: " +
                       std::to_string(paras.GetNMultiComplex()) +
                       " given but " + std::to_string(check_nMComplex) +
                       " expected."));
#endif

  /* Get parameters from ParameterList:
   * We use the same order of the parameters as was used during tree
   * construction */
  double mR = paras.GetDoubleParameter(0)->GetValue();
  double gA = paras.GetDoubleParameter(1)->GetValue();

  std::vector<std::complex<double>> gammaA =
      paras.GetMultiComplex(0)->GetValues();
  std::vector<double> barrierSqA = paras.GetMultiDouble(0)->GetValues();
  std::vector<std::complex<double>> gammaB =
      paras.GetMultiComplex(1)->GetValues();
  std::vector<double> barrierSqB = paras.GetMultiDouble(1)->GetValues();
  std::vector<std::complex<double>> gammaC =
      paras.GetMultiComplex(2)->GetValues();
  std::vector<double> barrierSqC = paras.GetMultiDouble(2)->GetValues();

  // invariant masses
  std::vector<double> mSq = paras.GetMultiDouble(3)->GetValues();

  std::vector<std::complex<double>> results(mSq.size(),
                                            std::complex<double>(0., 0.));
  // TODO make sure all vectors have to same size

  // calc function for each point
  for (unsigned int ele = 0; ele < mSq.size(); ele++) {
    try {
      /* Generally we need to add a factor q^{2J+1} to each channel term.
       * But since Flatte resonances are usually J=0 we neglect it here.*/
      results.at(ele) = AmpFlatteRes::dynamicalFunction(
          mSq.at(ele), mR, gA, gammaA.at(ele) * barrierSqA.at(ele),
          gammaB.at(ele) * barrierSqB.at(ele),
          gammaC.at(ele) * barrierSqC.at(ele));
    } catch (std::exception &ex) {
      LOG(error) << "FlatteStrategy::execute() | " << ex.what();
      throw(std::runtime_error("FlatteStrategy::execute() | "
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
