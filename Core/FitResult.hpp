//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//			Peter Weidenkaff
//-------------------------------------------------------------------------------
//! Optimizer Interface Base-Class.
/*! \class Optimizer
 * @file Optimizer.hpp
 * This class provides the interface to (external) optimization libraries or
 * routines. As it is pure virtual, one needs at least one implementation to
 * provide an optimizer for the analysis which varies free model-parameters. If
 * a new optimizer is derived from and fulfills this base-class, no change in
 * other modules are necessary to work with the new optimizer library or
 * routine.
 */

#ifndef _FITRESULT_HPP_
#define _FITRESULT_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/export.hpp>

#include "Core/AmpIntensity.hpp"
#include "Core/ParameterList.hpp"
#include "Core/TableFormater.hpp"
#include "Core/PhysConst.hpp"

namespace ComPWA {

class FitResult {
public:
  FitResult() : time(0), nSetsFractionError(0){};
  
  virtual ~FitResult(){};
  
  //! Set AmpIntensity.
  virtual void SetIntensity(std::shared_ptr<AmpIntensity> intens) {
    _intens = intens;
  }
  //! Set fraction list
  virtual void SetFractions(ParameterList ini) { fractionList = ini; }
  
  //! Set list of initial parameters
  virtual void SetInitialParameters(ParameterList iniPars) {
    initialParameters.DeepCopy(iniPars);
  }
  
  //! Set list of final fit parameters
  virtual void SetFinalParameters(ParameterList finPars);
  
  //! Set list of true parameters
  virtual void SetTrueParameters(ParameterList truePars) {
    trueParameters.DeepCopy(truePars);
  }
  
  //! Set value of likelihood with initial parameter
  virtual void SetInitialLH(double iniLH) {}
  
  //! Get list of initial parameters
  virtual ParameterList GetInitialParameters() { return initialParameters; }
  
  //! Get list of final fit parameters
  virtual ParameterList GetFinalParameters() { return finalParameters; }
  
  //! Get list of true parameters
  virtual ParameterList GetTrueParameters() { return trueParameters; }
  
  //! Set processing time for minimization
  virtual void SetTime(double t) { time = t; }
  
  //! Get processing time for minimization
  virtual double GetTime() { return time; }
  
  //! Get fit result (e.g. likelihood or chi2)
  virtual double GetResult() = 0;
  

  //! Table with fit parameters
  virtual void PrintFitParameters(TableFormater *tableResult);
  
  //! Table with fit fractions
  virtual void PrintFitFractions(TableFormater *tab);
  
  //! Table with fit fractions
  virtual void PrintFitFractions(TableFormater *tab,
                                 std::shared_ptr<AmpIntensity> amp,
                                 int nErrorSets = 0);
  
  //! Getter function for fractions list. Make sure that fractions are
  //! calculated beforehand.
  virtual ParameterList &GetFractions() { return fractionList; }

  //! Enable correct error estimation for fit fractions. Very time consuming!
  void SetUseCorrelatedErrors(int nSets = 200);

  //! Print fit result
  virtual void Print(std::string opt = "");

  virtual void WriteTeX(std::string filename){};
  virtual void WriteXML(std::string filename){};
  virtual void WriteText(std::string filename);
  virtual void WriteSimpleText(std::string filename);
  virtual operator double() const = 0;
  friend std::ostream &operator<<(std::ostream &out, FitResult &fitres) {
    out << fitres.GetResult();
    return out;
  };
  //! Any errors during minimization?
  virtual bool HasFailed() { return 0; };

protected:
  virtual double shiftAngle(double v);
  virtual void genOutput(std::ostream &out, std::string opt = "") = 0;
  virtual void genSimpleOutput(std::ostream &out);

  virtual double GetCorr(unsigned int n, unsigned int t) { return -9000; };
  
  //! Time for minimization
  double time;
  //! Initial list of parameters
  ParameterList initialParameters;
  //! Final list of parameters
  ParameterList finalParameters;
  //! True list of parameters
  ParameterList trueParameters;
  //! Fit amplitude (can't be serialized)
  std::shared_ptr<AmpIntensity> _intens;

  //! Number of parameter sets that are used to propagate the cov matrix through
  //! the normalization
  int nSetsFractionError;

  /** Calculate fit fractions.
   * Fractions are calculated using the formular:
   * \f[
   *  f_i = \frac{|c_i|^2 \int A_i A_i^*}{\int \sum c_l c_m^* A_l A_m}
   * \f]
   * The \f$c_i\f$ complex coefficienct of the amplitude and the denominatior is
   * the integral over
   * the whole amplitude.
   *
   * @param parList result with fit fractions for the single resonances
   * @param amp amplitude with the resonances which share the fit fractions
   * @param nSets Precise error calcucation using @param nSets Monte-Carlo
   * events
   */
  //! Calculate fit fractions and its errors.
  virtual void calcFractionError(ParameterList &parList,
                                 std::shared_ptr<AmpIntensity> amp,
                                 int nSets = 200) = 0;

  //! List with fit fractions and errors
  ParameterList fractionList;
  double sumFractions;
  double sumFractionsError;

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    ar &BOOST_SERIALIZATION_NVP(time);
    ar &BOOST_SERIALIZATION_NVP(initialParameters);
    ar &BOOST_SERIALIZATION_NVP(finalParameters);
    ar &BOOST_SERIALIZATION_NVP(trueParameters);
    ar &BOOST_SERIALIZATION_NVP(fractionList);
    ar &BOOST_SERIALIZATION_NVP(sumFractions);
    ar &BOOST_SERIALIZATION_NVP(sumFractionsError);
    ar &BOOST_SERIALIZATION_NVP(nSetsFractionError);
  }
};
BOOST_SERIALIZATION_ASSUME_ABSTRACT(FitResult);
} /* namespace ComPWA */

#endif
