// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Base class FitResult.
///

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
#include "Core/Properties.hpp"

namespace ComPWA {

class FitResult {
public:
  FitResult() : time(0), sumFractions(0.0), sumFractionsError(0.0) {};
  
  virtual ~FitResult(){};
  
  /// Set list of initial parameters
  virtual void SetInitialParameters(ParameterList iniPars) {
    initialParameters.DeepCopy(iniPars);
  }
  
  /// Get list of initial parameters
  virtual ParameterList GetInitialParameters() { return initialParameters; }
  
  /// Set list of final fit parameters
  virtual void SetFinalParameters(ParameterList &finPars);
  
  /// Get list of final fit parameters
  virtual ParameterList GetFinalParameters() { return finalParameters; }
  
  /// Set list of true parameters
  virtual void SetTrueParameters(ParameterList &truePars);
  
  /// Get list of true parameters
  virtual ParameterList GetTrueParameters() { return trueParameters; }
  
  /// Set processing time for minimization
  virtual void SetTime(double t) { time = t; }
  
  /// Get processing time for minimization
  virtual double GetTime() const { return time; }
  
  /// Get fit result (e.g. likelihood or chi2)
  virtual double GetResult() = 0;
  
  /// Set list with fit fractions
  virtual void SetFitFractions(ParameterList &list);
  
  /// Get list of fit fractions
  virtual ParameterList GetFitFractions() {
    return _fitFractions;
  }

  /// Table with fit parameters
  virtual void PrintFitParameters(TableFormater *tableResult);
  
  /// Table with fit fractions
  virtual void PrintFitFractions(TableFormater *tab);
  
  /// Print fit result
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
  /// Any errors during minimization?
  virtual bool HasFailed() { return 0; };

protected:
  virtual double shiftAngle(double v);
  virtual void genOutput(std::ostream &out, std::string opt = "") = 0;
  virtual void genSimpleOutput(std::ostream &out);

  virtual double GetCorr(unsigned int n, unsigned int t) { return -9000; };
  
  /// Time for minimization
  double time;
  
  /// Initial list of parameters
  ParameterList initialParameters;
  
  /// Final list of parameters
  ParameterList finalParameters;
  
  /// True list of parameters
  ParameterList trueParameters;
  
  ParameterList _fitFractions;
  
  /// List with fit fractions and errors
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
} // ns::ComPWA

#endif
