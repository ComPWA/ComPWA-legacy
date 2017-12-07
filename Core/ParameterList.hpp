// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// ParameterList class
///

#ifndef _PARAMETERLIST_HPP_
#define _PARAMETERLIST_HPP_

#include <iostream>
#include <string>
#include <memory>
#include <sstream>
#include <vector>
#include <map>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "Core/Parameter.hpp"
#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"
#include "Core/Value.hpp"

namespace ComPWA {

///
/// \class ParameterList
/// This class provides a list of parameters and values of different types.
///
class ParameterList {
public:
  ParameterList() {};
  /// Only shared_ptr are copied. Those still point to the same object.
  /// See DeepCopy(const ParameterList &in).
  ParameterList(const ParameterList &in) = default;

  /// Clear this parameter and deep-copy all parameters from \p in. Deep-copy
  /// means that for each parameter a new object is created (not only the
  /// shared_ptr is copied).
  void DeepCopy(const ParameterList &in){};

  virtual ~ParameterList(){};

  virtual std::size_t numParameters() const;

  virtual void addParameter(std::shared_ptr<Parameter> par);

  virtual void addParameters(std::vector<std::shared_ptr<Parameter>> pars);

  virtual std::size_t numValues() const;

  virtual void addValue(std::shared_ptr<Parameter> value);

  virtual void addValues(std::vector<std::shared_ptr<Parameter>> values);

  // Parameter
  virtual std::shared_ptr<DoubleParameter> doubleParameter(size_t i) {
    return DoubleParameters.at(i);
  };

  virtual std::vector<std::shared_ptr<DoubleParameter>> &doubleParameters() {
    return DoubleParameters;
  };

  // Value
  // Single sized values
  virtual std::shared_ptr<ComPWA::Value<int>> intValue(size_t i) {
    return IntValues.at(i);
  };

  virtual std::vector<std::shared_ptr<ComPWA::Value<int>>> &intValues() {
    return IntValues;
  };

  virtual std::shared_ptr<ComPWA::Value<double>> doubleValue(size_t i) {
    return DoubleValues.at(i);
  };

  virtual std::vector<std::shared_ptr<ComPWA::Value<double>>> &doubleValues() {
    return DoubleValues;
  };

  virtual std::shared_ptr<ComPWA::Value<std::complex<double>>>
  complexValue(size_t i) {
    return ComplexValues.at(i);
  };

  virtual std::vector<std::shared_ptr<ComPWA::Value<std::complex<double>>>> &
  complexValues() {
    return ComplexValues;
  };

  virtual std::shared_ptr<ComPWA::Value<std::vector<int>>> mIntValue(size_t i) {
    return MultiIntValues.at(i);
  };

  virtual std::vector<std::shared_ptr<ComPWA::Value<std::vector<int>>>> &
  mIntValues() {
    return MultiIntValues;
  };

  virtual std::shared_ptr<ComPWA::Value<std::vector<double>>>
  mDoubleValue(size_t i) {
    return MultiDoubleValues.at(i);
  };

  virtual std::vector<std::shared_ptr<ComPWA::Value<std::vector<double>>>> &
  mDoubleValues() {
    return MultiDoubleValues;
  };

  virtual std::shared_ptr<ComPWA::Value<std::vector<std::complex<double>>>>
  mComplexValue(size_t i) {
    return MultiComplexValues.at(i);
  };

  virtual std::vector<
      std::shared_ptr<ComPWA::Value<std::vector<std::complex<double>>>>> &
  mComplexValues() {
    return MultiComplexValues;
  };

protected:
  std::vector<std::shared_ptr<ComPWA::Value<int>>> IntValues;

  std::vector<std::shared_ptr<ComPWA::Value<double>>> DoubleValues;

  std::vector<std::shared_ptr<ComPWA::Value<std::complex<double>>>>
      ComplexValues;

  std::vector<std::shared_ptr<ComPWA::Value<std::vector<int>>>> MultiIntValues;

  std::vector<std::shared_ptr<ComPWA::Value<std::vector<double>>>>
      MultiDoubleValues;

  std::vector<std::shared_ptr<ComPWA::Value<std::vector<std::complex<double>>>>>
      MultiComplexValues;

  std::vector<std::shared_ptr<ComPWA::DoubleParameter>> DoubleParameters;

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    using namespace boost::serialization;
    // currently only DoubleParameters can be serialized
    ar &make_nvp("DoubleParameters", DoubleParameters);
  }
};

} // namespace ComPWA

/// Support for serialization of std::shared_ptr (and other types) is
/// added in boost 1.56 . For previous versions we have to implement
/// it ourself
#include <boost/version.hpp>
#if (BOOST_VERSION < 105600)
#include <boost/serialization/split_free.hpp>
#include <boost/unordered_map.hpp>
#include <typeinfo>

//--- Wrapper for std::shared_ptr<T> ------------------------------------------
namespace boost {
namespace serialization {

template <class Archive, class Type>
void save(Archive &archive, const std::shared_ptr<Type> &value,
          const unsigned int version) {
  Type *data = value.get();
  archive << make_nvp("shared_ptr", data);
}

template <class Archive, class Type>
void load(Archive &archive, std::shared_ptr<Type> &value,
          const unsigned int version) {
  Type *data;
  archive >> make_nvp("shared_ptr", data);
  //  archive >>data;

  typedef std::weak_ptr<Type> WeakPtr;
  static boost::unordered_map<void *, WeakPtr> hash;

  if (hash[data].expired()) {
    value = std::shared_ptr<Type>(data);
    hash[data] = value;
  } else
    value = hash[data].lock();
}

template <class Archive, class Type>
inline void serialize(Archive &archive, std::shared_ptr<Type> &value,
                      const unsigned int version) {
  split_free(archive, value, version);
}

} // ns:serialization
} // ns:boost
#endif // END serialization work-a-round

#endif
