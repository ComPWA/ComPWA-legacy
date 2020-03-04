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

#ifndef COMPWA_FUNCTIONS_HPP_
#define COMPWA_FUNCTIONS_HPP_

#include <complex>
#include <math.h>
#include <vector>

#include "Core/Exceptions.hpp"
#include "Core/FunctionTree/FitParameter.hpp"
#include "Core/FunctionTree/ParameterList.hpp"

namespace ComPWA {
namespace FunctionTree {
///
/// \class Strategy
/// Virtual base class for operations of FunctionTree nodes.
///
class Strategy {
public:
  Strategy(ParType out, std::string name) : checkType(out), Name(name){};

  virtual ~Strategy() = default;

  /// Return parameter type
  virtual const ParType OutType() const { return checkType; }

  /// Strategy execution
  virtual void execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out) = 0;

  std::string str() const { return Name; }

  friend std::ostream &operator<<(std::ostream &out,
                                  std::shared_ptr<Strategy> b) {
    return out << b->str();
  }

  friend std::ostream &operator<<(std::ostream &out, const Strategy &b) {
    return out << b.str();
  }

protected:
  ParType checkType;
  const std::string Name;
};

///
/// \class Inverse
/// Calculates the inverse of input double values and double parameters.
///
class Inverse : public Strategy {
public:
  Inverse(ParType out) : Strategy(out, "Inv"){};

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
  SquareRoot(ParType out) : Strategy(out, "Sqrt"){};

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
  AddAll(ParType out) : Strategy(out, "AddAll"){};

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
  MultAll(ParType out) : Strategy(out, "MultAll"){};

  virtual ~MultAll() {}

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class LogOf : public Strategy {
public:
  LogOf(ParType out) : Strategy(out, "LogOf"){};

  virtual ~LogOf(){};

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class Exp : public Strategy {
public:
  Exp(ParType out) : Strategy(out, "Exp"){};

  virtual ~Exp(){};

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class Pow : public Strategy {
  int power;

public:
  Pow(ParType out, int power_) : Strategy(out, "Pow"), power(power_){};

  virtual ~Pow(){};

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class Complexify : public Strategy {
public:
  Complexify(ParType out) : Strategy(out, "Complexify"){};

  virtual ~Complexify() {}

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class ComplexConjugate : public Strategy {
public:
  ComplexConjugate(ParType out) : Strategy(out, "ComplexConjugate"){};

  virtual ~ComplexConjugate() {}

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

class AbsSquare : public Strategy {
public:
  AbsSquare(ParType out) : Strategy(out, "AbsSquare"){};

  virtual ~AbsSquare() {}

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);
};

} // namespace FunctionTree
} // namespace ComPWA

#endif
