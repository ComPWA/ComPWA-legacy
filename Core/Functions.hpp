// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Functions to be used in FuntionTree.
/*! \class Strategy
 * \class AddAll
 * \class MultAll
 * \class PowerTwo
 * @file Functions.hpp
 * This file contains Functions implementing the Strategy interface so they
 * can be used inside a node of the FuntionTree to calculate the node-value.
 * In addition to the simple functions provided here, the interface can also
 * be used at other places to provide functions for the FunctionTree.
 */

#ifndef _FUNCTIONS_HPP_
#define _FUNCTIONS_HPP_

#include <vector>
#include <complex>
#include <math.h>

#include "Core/Exceptions.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"
#include "Core/AbsParameter.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

class Strategy {
public:
  //! Constructor
  Strategy(ParType in) : checkType(in){};
  virtual ~Strategy() {}

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type.
   * \sa make_str(), to_str()
   */
  friend std::ostream &operator<<(std::ostream &out,
                                  std::shared_ptr<Strategy> b) {
    return out << b->to_str();
  }

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type.
   * \sa make_str(), to_str()
   */
  friend std::ostream &operator<<(std::ostream &out, const Strategy &b) {
    return out << b.to_str();
  }

  //! Get ParType
  virtual const ParType OutType() const { return checkType; }

  //! Pure Virtual interface for streaming info about the strategy
  virtual const std::string to_str() const = 0;

  //! Pure Virtual interface for executing a strategy
  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out) = 0;

protected:
  ParType checkType;
};

class Inverse : public Strategy {
public:
  Inverse(ParType in) : Strategy(in){};
  virtual ~Inverse(){};

  virtual const std::string to_str() const { return "+"; }

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);
};

class SquareRoot : public Strategy {
public:
  SquareRoot(ParType in) : Strategy(in){};

  virtual ~SquareRoot() {}

  virtual const std::string to_str() const { return "+"; }

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);
};

class AddAll : public Strategy {
public:
  AddAll(ParType in) : Strategy(in){};
  virtual ~AddAll() {}

  virtual const std::string to_str() const { return "+"; }

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);
};

class MultAll : public Strategy {
public:
  MultAll(ParType in) : Strategy(in){};
  virtual ~MultAll() {}

  virtual const std::string to_str() const { return "*"; };

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);
};

class LogOf : public Strategy {
public:
  LogOf(ParType in) : Strategy(in){};
  virtual ~LogOf(){};

  virtual const std::string to_str() const { return "Log"; };

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);
};

class Real : public Strategy {
public:
  Real(ParType in) : Strategy(in){};
  virtual ~Real() {}

  virtual const std::string to_str() const { return "RealPart"; };

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out) {

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
      for (auto &complex_element :
           paras.GetMultiComplex(0)->GetValues()) {
        results.push_back(complex_element.real());
      }

      out = std::shared_ptr<AbsParameter>(
          new MultiDouble(out->GetName(), results));

      break;
    } // end multi complex

    case ParType::COMPLEX: {
      if (!(nC == 1)) {
        // TODO: exception wrong input
        return false;
      }

      out = std::shared_ptr<AbsParameter>(new DoubleParameter(
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
  Complexify(ParType in) : Strategy(in){};
  virtual ~Complexify() {}

  virtual const std::string to_str() const { return "MakeComplex"; };

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);
};

class ComplexConjugate : public Strategy {
public:
  ComplexConjugate(ParType in) : Strategy(in){};
  virtual ~ComplexConjugate() {}

  virtual const std::string to_str() const { return "ComplexConjugate"; };

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out) {
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
        out = std::shared_ptr<AbsParameter>(
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
      out = std::shared_ptr<AbsParameter>(new ComplexParameter(
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
  AbsSquare(ParType in) : Strategy(in){};
  virtual ~AbsSquare() {}

  virtual const std::string to_str() const { return "||^2"; };

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);
};

class Power : public Strategy {
public:
  Power(ParType in) : Strategy(in){};
  virtual ~Power() {}

  virtual const std::string to_str() const { return "^"; }

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);
};

} /* namespace ComPWA */

#endif
