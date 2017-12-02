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

template <typename T>
inline const char* typeName(void) { return "unknown"; }

template <>
inline const char* typeName<std::vector<std::complex<double>>>(void) { return "MComplex"; }

template <class T> class ParameterT : public ComPWA::Parameter {
public:
  ParameterT() : Parameter("") { std::cout << typeName<T>() << std::endl; }
  virtual std::string className() const {

    std::cout << typeName<T>() << std::endl;
  }
  /// A public function returning a string with parameter information
  virtual std::string to_str() const {};

  /// A public function returning a string with parameter value
  virtual std::string val_to_str() const {};
};

template class ParameterT<std::vector<std::complex<double>>>;

} // ns::ComPWA
#endif
