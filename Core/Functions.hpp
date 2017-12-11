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
  virtual void execute(ParameterList &paras,
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

///
/// \class Inverse
/// Calculates the inverse of input double values and double parameters.
///
class Inverse : public Strategy {
public:
  Inverse(ParType in) : Strategy(in, "Inv"){};

  virtual ~Inverse(){};

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

///
/// \class SquareRoot
/// Calculates the square root of input double values and double parameters.
/// Complex parameters are currently not supported.
///
class SquareRoot : public Strategy {
public:
  SquareRoot(ParType in) : Strategy(in, "Sqrt"){};

  virtual ~SquareRoot() {}

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

///
/// \class AddAll
/// Calculates the square root of input double values and double parameters.
/// Complex parameters are currently not supported.
///
class AddAll : public Strategy {
public:
  AddAll(ParType in) : Strategy(in, "Add"){};

  virtual ~AddAll() {}

/// Add all values. Depending on the output type the summation is a little
/// different:
///   - ParType::MCOMPLEX: all single complex, integer and double values are
///     added to a std::complex<double>. Each multi value is added element by
///     element and the previous result from the single values is added to
///     each element.
///   - ParType::MDOUBLE: same ad MCOMPLEX except that complex 
  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class MultAll : public Strategy {
public:
  MultAll(ParType in) : Strategy(in, "Mult"){};

  virtual ~MultAll() {}

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class LogOf : public Strategy {
public:
  LogOf(ParType in) : Strategy(in, "LogOf"){};

  virtual ~LogOf(){};

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class Complexify : public Strategy {
public:
  Complexify(ParType in) : Strategy(in, "Complexify"){};

  virtual ~Complexify() {}

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class ComplexConjugate : public Strategy {
public:
  ComplexConjugate(ParType in) : Strategy(in, "ComplexConjugate"){};

  virtual ~ComplexConjugate() {}
  
  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class AbsSquare : public Strategy {
public:
  AbsSquare(ParType in) : Strategy(in, "Abs"){};

  virtual ~AbsSquare() {}

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

} // ns::ComPWA

#endif
