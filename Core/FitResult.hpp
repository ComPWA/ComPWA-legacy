// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

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
  FitResult() : time(0) {};
  
  virtual ~FitResult(){};
  
  //! Set AmpIntensity.
  virtual void SetIntensity(std::shared_ptr<AmpIntensity> intens) {
    _intens = intens;
  }
  
  //! Set list of initial parameters
  virtual void SetInitialParameters(ParameterList iniPars) {
    initialParameters.DeepCopy(iniPars);
  }
  
  //! Get list of initial parameters
  virtual ParameterList GetInitialParameters() { return initialParameters; }
  
  //! Set list of final fit parameters
  virtual void SetFinalParameters(ParameterList finPars);
  
  //! Get list of final fit parameters
  virtual ParameterList GetFinalParameters() { return finalParameters; }
  
  //! Set list of true parameters
  virtual void SetTrueParameters(ParameterList truePars) {
    trueParameters.DeepCopy(truePars);
  }
  
  //! Get list of true parameters
  virtual ParameterList GetTrueParameters() { return trueParameters; }
  
  //! Set processing time for minimization
  virtual void SetTime(double t) { time = t; }
  
  //! Get processing time for minimization
  virtual double GetTime() const { return time; }
  
  //! Get fit result (e.g. likelihood or chi2)
  virtual double GetResult() = 0;
  
  //! Set list with fit fractions
  virtual void SetFitFractions(ParameterList list) {
    _fitFractions.DeepCopy(list);
  }
  
  //! Get list of fit fractions
  virtual ParameterList GetFitFractions() {
    return _fitFractions;
  }

  //! Table with fit parameters
  virtual void PrintFitParameters(TableFormater *tableResult);
  
  //! Table with fit fractions
  virtual void PrintFitFractions(TableFormater *tab);
  
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

  ParameterList _fitFractions;
  
  //! List with fit fractions and errors
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
    ar &BOOST_SERIALIZATION_NVP(_fitFractions);
    ar &BOOST_SERIALIZATION_NVP(sumFractions);
    ar &BOOST_SERIALIZATION_NVP(sumFractionsError);
  }
};
BOOST_SERIALIZATION_ASSUME_ABSTRACT(FitResult);
} /* namespace ComPWA */

#endif
