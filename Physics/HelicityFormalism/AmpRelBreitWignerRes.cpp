//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - correct nominator, using dataPoint for data
// handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include <stdlib.h>

#include "Physics/HelicityFormalism/AmpRelBreitWignerRes.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

AmpRelBreitWignerRes::~AmpRelBreitWignerRes() {}

std::string AmpRelBreitWignerRes::to_str() const {
  std::string dynAmp; //= AmpAbsDynamicalFunction::to_str();
  std::stringstream str;
  str << _width->to_str() << std::endl;

  return dynAmp + str.str();
}

void AmpRelBreitWignerRes::CheckModified() const {
  AbstractDynamicalFunction::CheckModified();
  if (_width->GetValue() != _current_width) {
    SetModified();
    const_cast<double &>(_current_width) = _width->GetValue();
  }
  return;
}

double AmpRelBreitWignerRes::GetNormalization() const {
  CheckModified();
  if( GetModified() ) {
    const_cast<double &>(_current_integral) = Integral();
    SetModified(false);
  }
  return (1/_current_integral);
}

std::complex<double> AmpRelBreitWignerRes::Evaluate(const dataPoint &point,
                                                    int pos) const {
  double mSq = point.GetValue(pos);
  std::complex<double> result;
  try {
    result = dynamicalFunction(mSq, _mass->GetValue(), _massA, _massB,
                               _width->GetValue(), (double)_spin,
                               _mesonRadius->GetValue(), _ffType);
  } catch (std::exception &ex) {
    LOG(error) << "AmpRelBreitWignerRes::EvaluateAmp() | "
                  "Dynamical function can not be evalutated: "
               << ex.what();
    throw;
  }
  return result;
}

std::complex<double> AmpRelBreitWignerRes::dynamicalFunction(
    double mSq, double mR, double ma, double mb, double width, unsigned int J,
    double mesonRadius, formFactorType ffType) {
  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  double barrier =
      HelicityKinematics::FormFactor(sqrtS, ma, mb, J, mesonRadius, ffType) /
      HelicityKinematics::FormFactor(mR, ma, mb, J, mesonRadius, ffType);

  std::complex<double> qTerm = std::pow((Kinematics::phspFactor(sqrtS, ma, mb) /
                                         Kinematics::phspFactor(mR, ma, mb)) *
                                            mR / sqrtS,
                                        (2 * J + 1));

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

#ifndef NDEBUG
  if (std::isnan(result.real()) || std::isnan(result.imag())) {
    std::cout << "AmpRelBreitWignerRes::dynamicalFunction() | " << barrier
              << " " << mR << " " << mSq << " " << ma << " " << mb << std::endl;
    return 0;
  }
#endif

  return result;
}

std::shared_ptr<FunctionTree>
AmpRelBreitWignerRes::GetTree(ParameterList &sample,
                              ParameterList &toySample, int pos, std::string suffix) {
  return std::shared_ptr<FunctionTree>();
  
//    auto kin = dynamic_cast<DalitzKinematics *>(Kinematics::instance());
//    //	auto var1_limit = kin->GetMinMax( GetVarIdA() );
//    //	auto var2_limit = kin->GetMinMax( GetVarIdB() );
//    //	double phspVol = (var1_limit.second-var1_limit.first)
//    //			*(var2_limit.second-var2_limit.first);
//    double phspVol = kin->GetPhspVolume();
//    //	double phspVol = 1;
//  
//    int sampleSize = sample.GetMultiDouble(0)->GetNValues();
//    int toySampleSize = toySample.GetMultiDouble(0)->GetNValues();
//  
//    LOG(info) << "AmpRelBreitWignerRes::setupBasicTree() | " << _name
//              << " nEvents=" << sampleSize << " nPhspEvents=" <<
//              toySampleSize;
//  
//    //------------Setup Tree---------------------
//    std::shared_ptr<FunctionTree> newTree(new FunctionTree());
//  
//    //----Strategies needed
//    std::shared_ptr<MultAll> mmultStrat(new MultAll(ParType::MCOMPLEX));
//    std::shared_ptr<AbsSquare> msqStrat(new AbsSquare(ParType::MDOUBLE));
//    std::shared_ptr<MultAll> multDStrat(new MultAll(ParType::DOUBLE));
//    std::shared_ptr<AddAll> addStrat(new AddAll(ParType::DOUBLE));
//    std::shared_ptr<Complexify> complStrat(new Complexify(ParType::COMPLEX));
//    std::shared_ptr<Inverse> invStrat(new Inverse(ParType::DOUBLE));
//    std::shared_ptr<SquareRoot> sqRootStrat(new SquareRoot(ParType::DOUBLE));
//  
//    //----Add Nodes
//    std::shared_ptr<BreitWignerStrategy> rbwStrat(new
//    BreitWignerStrategy(_name));
//    std::shared_ptr<WignerDStrategy> angdStrat(new WignerDStrategy(_name));
//  
//    // Reso=BW*C*AD*N
//    newTree->createHead("Reso_" + _name, mmultStrat, sampleSize);
//  
//    newTree->createNode("PreFactor_" + _name, complStrat, "Reso_" + _name);
//    newTree->createLeaf("IntensPre_" + _name, std::abs(_preFactor),
//                        "PreFactor_" + _name);
//    newTree->createLeaf("PhasePre_" + _name, std::arg(_preFactor),
//                        "PreFactor_" + _name);
//  
//    newTree->createNode("C_" + _name, complStrat, "Reso_" + _name); // m0c
//    newTree->createLeaf("Intens_" + _name, _mag, "C_" + _name);     // r
//    newTree->createLeaf("Phase_" + _name, _phase, "C_" + _name);    // phi
//    // Angular distribution
//    if ((double)_spin)
//      newTree->insertTree(_wignerD.SetupTree(sample, suffix), "Reso_" +
//      _name);
//  
//    // Breit-Wigner
//    newTree->createNode("RelBW_" + _name, rbwStrat, "Reso_" + _name,
//    sampleSize);
//    newTree->createLeaf("mass", _mass, "RelBW_" + _name);         // m0
//    newTree->createLeaf("width", _width, "RelBW_" + _name);       // resWidth
//    newTree->createLeaf("spin", (double)_spin, "RelBW_" + _name); // spin
//    newTree->createLeaf("mesonRadius", _mesonRadius, "RelBW_" + _name); // d
//    newTree->createLeaf("formFactorType", _ffType, "RelBW_" + _name);   // d
//    newTree->createLeaf("ma", _mass1, "RelBW_" + _name);                // ma
//    newTree->createLeaf("mb", _mass2, "RelBW_" + _name);                // mb
//    newTree->createLeaf("sample", sample.GetMultiDouble(_subSys),
//                        "RelBW_" + _name); // mc
//  
//    // Normalization
//    if (_normStyle == normStyle::none) {
//      newTree->createLeaf("N_" + _name, 1., "Reso_" + _name);
//    } else {
//      // Normalization parameter for dynamical amplitude
//      newTree->createNode("N_" + _name, sqRootStrat,
//                          "Reso_" + _name); // N = sqrt(NSq)
//      newTree->createNode(
//          "NSq_" + _name, multDStrat,
//          "N_" + _name); // NSq = (PhspVolume/N_phspMC * Sum(|A|^2))^-1
//      newTree->createLeaf("PhspSize_" + _name, toySampleSize,
//                          "NSq_" + _name); // N_phspMC
//      newTree->createLeaf("PhspVolume_" + _name, 1 / phspVol,
//                          "NSq_" + _name); // 1/PhspVolume
//      newTree->createNode("InvSum_" + _name, invStrat,
//                          "NSq_" + _name); // 1/Sum(|A|^2)
//      newTree->createNode("Sum_" + _name, addStrat,
//                          "InvSum_" + _name); // Sum(|A|^2)
//      newTree->createNode("AbsVal_" + _name, msqStrat, "Sum_" + _name);
//      //|A_i|^2
//  
//      newTree->createNode("NormReso_" + _name, mmultStrat, "AbsVal_" + _name,
//                          toySampleSize);
//  
//      // Angular distribution (Normalization)
//      if ((double)_spin)
//        newTree->insertTree(_wignerD.SetupTree(toySample, suffix + "_norm"),
//                            "NormReso_" + _name);
//      // Breit-Wigner (Normalization)
//      newTree->createNode("NormBW_" + _name, rbwStrat, "NormReso_" + _name,
//                          toySampleSize);                            // BW
//      newTree->createLeaf("mass", _mass, "NormBW_" + _name);         // m0
//      newTree->createLeaf("width", _width, "NormBW_" + _name);       //
//      resWidth
//      newTree->createLeaf("spin", (double)_spin, "NormBW_" + _name); // spin
//      newTree->createLeaf("mesonRadius", _mesonRadius, "NormBW_" + _name); //
//      d
//      newTree->createLeaf("formFactorType", _ffType, "NormBW_" + _name);   //
//      d
//      newTree->createLeaf("ma", _mass1, "NormBW_" + _name);                //
//      ma
//      newTree->createLeaf("mb", _mass2, "NormBW_" + _name);                //
//      mb
//      newTree->createLeaf("phspSample", toySample.GetMultiDouble(_subSys),
//                          "NormBW_" + _name);
//    }
//    return newTree;
}

//bool BreitWignerStrategy::execute(ParameterList &paras,
//                                  std::shared_ptr<AbsParameter> &out) {
//#ifndef NDEBUG
//  // Check parameter type
//  if (checkType != out->type())
//    throw(WrongParType("BreitWignerStrat::execute() | "
//                       "Output parameter is of type " +
//                       std::string(ParNames[out->type()]) +
//                       " and conflicts with expected type " +
//                       std::string(ParNames[checkType])));
//
//  // How many parameters do we expect?
//  int check_nBool = 0;
//  int check_nInt = 0;
//  int check_nComplex = 0;
//  int check_nDouble = 7;
//  int check_nMDouble = 1;
//  int check_nMComplex = 0;
//
//  // Check size of parameter list
//  if (paras.GetNBool() != check_nBool)
//    throw(BadParameter("BreitWignerStrat::execute() | "
//                       "Number of BoolParameters does not match: " +
//                       std::to_string(paras.GetNBool()) + " given but " +
//                       std::to_string(check_nBool) + " expected."));
//  if (paras.GetNInteger() != check_nInt)
//    throw(BadParameter("BreitWignerStrat::execute() | "
//                       "Number of IntParameters does not match: " +
//                       std::to_string(paras.GetNInteger()) + " given but " +
//                       std::to_string(check_nInt) + " expected."));
//  if (paras.GetNDouble() != check_nDouble)
//    throw(BadParameter("BreitWignerStrat::execute() | "
//                       "Number of DoubleParameters does not match: " +
//                       std::to_string(paras.GetNDouble()) + " given but " +
//                       std::to_string(check_nDouble) + " expected."));
//  if (paras.GetNComplex() != check_nComplex)
//    throw(BadParameter("BreitWignerStrat::execute() | "
//                       "Number of ComplexParameters does not match: " +
//                       std::to_string(paras.GetNComplex()) + " given but " +
//                       std::to_string(check_nComplex) + " expected."));
//  if (paras.GetNMultiDouble() != check_nMDouble)
//    throw(BadParameter("BreitWignerStrat::execute() | "
//                       "Number of MultiDoubles does not match: " +
//                       std::to_string(paras.GetNMultiDouble()) + " given but " +
//                       std::to_string(check_nMDouble) + " expected."));
//  if (paras.GetNMultiComplex() != check_nMComplex)
//    throw(BadParameter("BreitWignerStrat::execute() | "
//                       "Number of MultiComplexes does not match: " +
//                       std::to_string(paras.GetNMultiComplex()) +
//                       " given but " + std::to_string(check_nMComplex) +
//                       " expected."));
//#endif
//
//  /** Get parameters from ParameterList:
//   * We use the same order of the parameters as was used during tree
//   * construction
//   */
//  double m0 = paras.GetDoubleParameter(0)->GetValue();
//  double Gamma0 = paras.GetDoubleParameter(1)->GetValue();
//  unsigned int spin = (unsigned int)paras.GetDoubleParameter(2)->GetValue();
//  double d = paras.GetDoubleParameter(3)->GetValue();
//  formFactorType ffType =
//      formFactorType(paras.GetDoubleParameter(4)->GetValue());
//  double ma = paras.GetDoubleParameter(5)->GetValue();
//  double mb = paras.GetDoubleParameter(6)->GetValue();
//
//  std::vector<double> mp = paras.GetMultiDouble(0)->GetValues();
//
//  std::vector<std::complex<double>> results(mp.size(),
//                                            std::complex<double>(0., 0.));
//
//  // calc function for each point
//  for (unsigned int ele = 0; ele < mp.size(); ele++) {
//    try {
//      results.at(ele) = AmpRelBreitWignerRes::dynamicalFunction(
//          mp.at(ele), m0, ma, mb, Gamma0, spin, d, ffType);
//    } catch (std::exception &ex) {
//      LOG(error) << "BreitWignerStrategy::execute() | " << ex.what();
//      throw(std::runtime_error("BreitWignerStrategy::execute() | "
//                               "Evaluation of dynamic function failed!"));
//    }
//  }
//  out =
//      std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(), results));
//  return true;
//}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
