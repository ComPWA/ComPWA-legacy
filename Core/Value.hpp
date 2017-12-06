// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Template implementation of Parameter.
///

#ifndef ParameterT_hpp
#define ParameterT_hpp

#include <stdio.h>
#include "Core/Parameter.hpp"
namespace ComPWA {


template <class T> class Value : public ComPWA::Parameter {
public:
  Value() : Parameter("") { Type = typeName<T>(); }

  virtual T value() const { return Val; }

  virtual void setValue(T inVal) { Val = inVal; };

  operator T() const { return Val; };

  bool operator==(const BoolParameter otherPar) const;
  /// A public function returning a string with parameter information
  virtual std::string to_str() const { return ""; };

  /// A public function returning a string with parameter value
  virtual std::string val_to_str() const { return ""; };

protected:
  virtual std::string className() const {
    return std::string(ParNames[(int)Type]);
  }

  T Val;
};


} // ns::ComPWA
#endif
