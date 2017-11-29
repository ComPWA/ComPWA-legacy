// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file Functions.hpp
/// This file contains Functions implementing the Strategy interface so they
/// can be used inside a node of the FuntionTree to calculate the node-value.
/// In addition to the simple functions provided here, the interface can also
/// be used at other places to provide functions for the FunctionTree.
///

#ifndef _FUNCTIONS_HPP_
#define _FUNCTIONS_HPP_

#include <vector>
#include <complex>
#include <math.h>

#include "Core/Exceptions.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

///
/// \class Strategy
/// Virtual base class for operations of FunctionTree nodes.
///
class Strategy {
public:
  Strategy(ParType in, std::string op = "") : checkType(in), Op(op){};

  virtual ~Strategy() {}

  /// Return parameter type
  virtual const ParType OutType() const { return checkType; }

  /// Strategy execution
  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out) = 0;

  std::string str() const { return Op; }
  
  friend std::ostream &operator<<(std::ostream &out,
                                  std::shared_ptr<Strategy> b) {
    return out << b->str();
  }

  friend std::ostream &operator<<(std::ostream &out, const Strategy &b) {
    return out << b.str();
  }
  
protected:
  ParType checkType;
  const std::string Op;

};

class Inverse : public Strategy {
public:
  Inverse(ParType in) : Strategy(in, "Inv"){};

  virtual ~Inverse(){};

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out);
};

class SquareRoot : public Strategy {
public:
  SquareRoot(ParType in) : Strategy(in, "Sqrt"){};

  virtual ~SquareRoot() {}

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out);
};

class AddAll : public Strategy {
public:
  AddAll(ParType in) : Strategy(in, "Add"){};

  virtual ~AddAll() {}

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out);
};

class MultAll : public Strategy {
public:
  MultAll(ParType in) : Strategy(in, "Mult"){};

  virtual ~MultAll() {}

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out);
};

class LogOf : public Strategy {
public:
  LogOf(ParType in) : Strategy(in, "LogOf"){};

  virtual ~LogOf(){};

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out);
};

class Real : public Strategy {
public:
  Real(ParType in) : Strategy(in, "Real"){};

  virtual ~Real() {}

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out) {

    unsigned int nMC = paras.GetNMultiComplex();
    unsigned int nC = paras.GetNComplex();

    switch (checkType) {

    case ParType::MCOMPLEX: {
      if (!(nMC == 1)) {
        // TODO: exception wrong input
        return false;
      }
      unsigned int nElements = paras.GetMultiComplex(0)->GetNValues();
      // fill MultiDouble parameter
      std::vector<double> results;
      results.reserve(nElements);
      for (auto const &complex_element :
           paras.GetMultiComplex(0)->GetValues()) {
        results.push_back(complex_element.real());
      }

      out = std::shared_ptr<Parameter>(
          new MultiDouble(out->GetName(), results));

      break;
    } // end multi complex

    case ParType::COMPLEX: {
      if (!(nC == 1)) {
        // TODO: exception wrong input
        return false;
      }

      out = std::shared_ptr<Parameter>(new DoubleParameter(
          out->GetName(), paras.GetComplexParameter(0)->GetValue().real()));

      break;
    } // end double

    default: {
      // TODO: exception output partype wrong
      return false;
    }

    } // end switch

    return true;
  };
};

class Complexify : public Strategy {
public:
  Complexify(ParType in) : Strategy(in, "Complexify"){};

  virtual ~Complexify() {}

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out);
};

class ComplexConjugate : public Strategy {
public:
  ComplexConjugate(ParType in) : Strategy(in, "ComplexConjugate"){};

  virtual ~ComplexConjugate() {}

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out) {
    if (checkType != out->type())
      return false;

    unsigned int nMC = paras.GetNMultiComplex();
    unsigned int nC = paras.GetNComplex();

    if (nMC + nC == 0) {
      // TODO: exception no input
      return false;
    }

    switch (checkType) {
    case ParType::MCOMPLEX: {
      // output complex: input must be one multicomplex
      if (!(nMC == 1)) {
        // TODO: exception wrong input
        return false;
      }

      for (auto const &multi_complex : paras.GetMultiComplexs()) {
        unsigned int nElements = multi_complex->GetNValues();
        // fill MultiDouble parameter
        std::vector<std::complex<double>> results;
        results.reserve(nElements);
        for (unsigned int ele = 0; ele < nElements; ele++) {
          results.push_back(std::conj(multi_complex->GetValue(ele)));
        }
        out = std::shared_ptr<Parameter>(
            new MultiComplex(out->GetName(), results));
      }
      break;
    } // end multi complex

    case ParType::COMPLEX: {
      // output complex: input must be a complex
      if (!(nC == 1)) {
        // TODO: exception wrong input
        return false;
      }
      out = std::shared_ptr<Parameter>(new ComplexParameter(
          out->GetName(), std::conj(paras.GetComplexParameter(0)->GetValue())));
      break;
    } // end double

    default: {
      // TODO: exception output partype wrong
      return false;
    }

    } // end switch

    return true;
  };
};

class AbsSquare : public Strategy {
public:
  AbsSquare(ParType in) : Strategy(in, "Abs"){};

  virtual ~AbsSquare() {}

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out);
};

class Power : public Strategy {
public:
  Power(ParType in) : Strategy(in, "PowerTo"){};

  virtual ~Power() {}

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out);
};

} // ns::ComPWA

#endif
