// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Template implementation of Parameter for simple values.
///

#ifndef ParameterT_hpp
#define ParameterT_hpp

#include "FitParameter.hpp"
#include <iterator>
namespace ComPWA {
namespace FunctionTree {

template <class T>
std::ostream &operator<<(std::ostream &stream, const std::vector<T> &values) {
  auto n = values.size();
  if (n > 5)
    n = 5; // Print first 5 elements
  auto first = values.begin();
  auto last = values.begin();
  std::advance(last, n);
  std::copy(first, last, std::ostream_iterator<T>(stream, ", "));
  if (values.size() > n)
    stream << "...";

  return stream;
}

template <class T> class Value : public Parameter {
public:
  Value(std::string name = "") : Parameter(name), Val(0) {
    Type = typeName<T>();
  }

  Value(T val) : Parameter(""), Val(val) { Type = typeName<T>(); }

  Value(std::string name, T val) : Parameter(name), Val(val) {
    Type = typeName<T>();
  }

  virtual T value() const { return Val; }

  /// Reference on the value. In case of T = std::vector<T2> a reference to the
  /// vector is returned.
  virtual T &values() { return Val; }

  virtual void setValue(T inVal) {
    Val = inVal;
    notify();
  };

  /// Conversion operator for internal type
  operator T() const { return Val; };

  const T &operator()() const { return Val; };

  T &operator()() { return Val; };

  /// A public function returning a string with parameter information
  virtual std::string to_str() const {
    return Name + " [" + ParNames[Type] + "]";
  };

  /// A public function returning a string with parameter value
  virtual std::string val_to_str() const {
    std::stringstream stream;
    stream << Val;
    return stream.str();
  };

protected:
  virtual std::string className() const {
    return std::string(ParNames[(int)Type]);
  }

  T Val;
};

inline std::shared_ptr<Parameter> ValueFactory(ParType t,
                                               std::string name = "") {
  std::shared_ptr<Parameter> p;
  switch (t) {
  case ParType::MCOMPLEX: {
    p = std::make_shared<Value<std::vector<std::complex<double>>>>(name);
    break;
  }
  case ParType::MDOUBLE: {
    p = std::make_shared<Value<std::vector<double>>>(name);
    break;
  }
  case ParType::MINTEGER: {
    p = std::make_shared<Value<std::vector<int>>>(name);
    break;
  }
  case ParType::COMPLEX: {
    p = std::make_shared<Value<std::complex<double>>>(name);
    break;
  }
  case ParType::DOUBLE: {
    p = std::make_shared<Value<double>>(name);
    break;
  }
  case ParType::INTEGER: {
    p = std::make_shared<Value<int>>(name);
    break;
  }
  default: {
    throw BadParameter("ValueFactory() | Parameter type " + std::to_string(t) +
                       " unknown!");
  }
  }
  return p;
}

inline std::shared_ptr<Value<std::vector<std::complex<double>>>>
MComplex(std::string name, size_t s,
         std::complex<double> el = std::complex<double>(0., 0.)) {

  return std::make_shared<Value<std::vector<std::complex<double>>>>(
      name, std::vector<std::complex<double>>(s, el));
}

inline std::shared_ptr<Value<std::vector<std::complex<double>>>>
MComplex(std::string name, std::vector<std::complex<double>> v) {

  return std::make_shared<Value<std::vector<std::complex<double>>>>(name, v);
}

inline std::shared_ptr<Value<std::vector<double>>>
MDouble(std::string name, size_t s, double el = 0.) {

  return std::make_shared<Value<std::vector<double>>>(
      name, std::vector<double>(s, el));
}

inline std::shared_ptr<Value<std::vector<double>>>
MDouble(std::string name, std::vector<double> v) {

  return std::make_shared<Value<std::vector<double>>>(name, v);
}

inline std::shared_ptr<Value<std::vector<int>>>
MInteger(std::string name, size_t s, int el = 0.) {

  return std::make_shared<Value<std::vector<int>>>(name,
                                                   std::vector<int>(s, el));
}

inline std::shared_ptr<Value<std::vector<int>>> MInteger(std::string name,
                                                         std::vector<int> v) {

  return std::make_shared<Value<std::vector<int>>>(name, v);
}

} // namespace FunctionTree
} // namespace ComPWA
#endif
