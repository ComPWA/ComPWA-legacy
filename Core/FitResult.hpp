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

#include "Core/ParameterList.hpp"
#include "Core/TableFormater.hpp"
#include "Core/Properties.hpp"

namespace ComPWA {

class FitResult {
public:
  FitResult() : Time(0), SumFractions(0.0), SumFractionsError(0.0) {};
  
  virtual ~FitResult(){};
  
  /// Set list of initial parameters
  virtual void setInitialParameters(ParameterList iniPars) {
    InitialParameters.DeepCopy(iniPars);
  }
  
  /// Get list of initial parameters
  virtual ParameterList initialParameters() { return InitialParameters; }
  
  /// Set list of final fit parameters
  virtual void setFinalParameters(ParameterList &finPars);
  
  /// Get list of final fit parameters
  virtual ParameterList finalParameters() { return FinalParameters; }
  
  /// Set list of true parameters
  virtual void setTrueParameters(ParameterList &truePars);
  
  /// Get list of true parameters
  virtual ParameterList trueParameters() { return TrueParameters; }
  
  /// Set processing time for minimization
  virtual void setTime(double t) { Time = t; }
  
  /// Get processing time for minimization
  virtual double time() const { return Time; }
  
  /// Get fit result (e.g. likelihood or chi2)
  virtual double result() = 0;
  
  /// Set list with fit fractions
  virtual void setFitFractions(ParameterList &list);
  
  /// Get list of fit fractions
  virtual ParameterList fitFractions() {
    return FitFractions;
  }

  /// Table with fit parameters
  virtual void printFitParameters(TableFormater *tableResult);
  
  /// Table with fit fractions
  virtual void printFitFractions(TableFormater *tab);
  
  /// Print fit result
  virtual void print(std::string opt = "");

  virtual void writeTeX(std::string filename){};
  virtual void writeXML(std::string filename){};
  virtual void writeText(std::string filename);
  virtual void writeSimpleText(std::string filename);
  virtual operator double() const = 0;
  friend std::ostream &operator<<(std::ostream &out, FitResult &fitres) {
    out << fitres.result();
    return out;
  };
  
  /// Any errors during minimization?
  virtual bool hasFailed() { return 0; };

protected:
  virtual double shiftAngle(double v);
  virtual void genOutput(std::ostream &out, std::string opt = "") = 0;
  virtual void genSimpleOutput(std::ostream &out);

  virtual double GetCorr(unsigned int n, unsigned int t) { return -9000; };
  
  /// Time for minimization
  double Time;
  
  /// Initial list of parameters
  ParameterList InitialParameters;
  
  /// Final list of parameters
  ParameterList FinalParameters;
  
  /// True list of parameters
  ParameterList TrueParameters;
  
  ParameterList FitFractions;
  
  /// List with fit fractions and errors
  double SumFractions;
  double SumFractionsError;

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    ar &BOOST_SERIALIZATION_NVP(Time);
    ar &BOOST_SERIALIZATION_NVP(InitialParameters);
    ar &BOOST_SERIALIZATION_NVP(FinalParameters);
    ar &BOOST_SERIALIZATION_NVP(TrueParameters);
    ar &BOOST_SERIALIZATION_NVP(FitFractions);
    ar &BOOST_SERIALIZATION_NVP(SumFractions);
    ar &BOOST_SERIALIZATION_NVP(SumFractionsError);
  }
};
BOOST_SERIALIZATION_ASSUME_ABSTRACT(FitResult);
} // ns::ComPWA

#endif
