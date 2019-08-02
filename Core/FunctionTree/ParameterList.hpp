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
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

#include "Core/Exceptions.hpp"
#include "Core/FunctionTree/FitParameter.hpp"
#include "Core/FunctionTree/Value.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {
namespace Data {
struct DataSet;
}
namespace FunctionTree {
///
/// \class ParameterList
/// This class provides a list of parameters and values of different types.
///
class ParameterList {
public:
  ParameterList(){};
  ParameterList(const ComPWA::Data::DataSet &DataSample);
  /// Only shared_ptr are copied. Those still point to the same object.
  /// See DeepCopy(const ParameterList &in).
  ParameterList(const ParameterList &in) = default;

  /// Clear this parameter and deep-copy all parameters from \p in. Deep-copy
  /// means that for each parameter a new object is created (not only the
  /// shared_ptr is copied).
  void DeepCopy(const ParameterList &in);

  virtual ~ParameterList(){};

  virtual std::size_t numParameters() const;

  std::shared_ptr<FitParameter>
  addUniqueParameter(std::shared_ptr<FitParameter> par);

  virtual void addParameter(std::shared_ptr<Parameter> par);

  virtual void addParameter(std::shared_ptr<FitParameter> par);

  virtual void addParameters(std::vector<std::shared_ptr<Parameter>> pars);

  virtual std::size_t numValues() const;

  virtual void addValue(std::shared_ptr<Parameter> value);

  virtual void addValues(std::vector<std::shared_ptr<Parameter>> values);

  // Parameter
  virtual std::shared_ptr<FitParameter> doubleParameter(size_t i) const {
    return FitParameters.at(i);
  };

  virtual std::vector<std::shared_ptr<FitParameter>> &doubleParameters() {
    return FitParameters;
  };

  virtual const std::vector<std::shared_ptr<FitParameter>> &
  doubleParameters() const {
    return FitParameters;
  };

  // Value
  // Single sized values
  virtual std::shared_ptr<Value<int>> intValue(size_t i) {
    return IntValues.at(i);
  };

  virtual std::vector<std::shared_ptr<Value<int>>> &intValues() {
    return IntValues;
  };

  virtual const std::vector<std::shared_ptr<Value<int>>> &intValues() const {
    return IntValues;
  };

  virtual std::shared_ptr<Value<double>> doubleValue(size_t i) const {
    return DoubleValues.at(i);
  };

  virtual std::vector<std::shared_ptr<Value<double>>> &doubleValues() {
    return DoubleValues;
  };

  virtual const std::vector<std::shared_ptr<Value<double>>> &
  doubleValues() const {
    return DoubleValues;
  };

  virtual std::shared_ptr<Value<std::complex<double>>>
  complexValue(size_t i) const {
    return ComplexValues.at(i);
  };

  virtual std::vector<std::shared_ptr<Value<std::complex<double>>>> &
  complexValues() {
    return ComplexValues;
  };

  virtual const std::vector<std::shared_ptr<Value<std::complex<double>>>> &
  complexValues() const {
    return ComplexValues;
  };

  virtual std::shared_ptr<Value<std::vector<int>>> mIntValue(size_t i) const {
    return MultiIntValues.at(i);
  };

  virtual std::vector<std::shared_ptr<Value<std::vector<int>>>> &mIntValues() {
    return MultiIntValues;
  };

  virtual const std::vector<std::shared_ptr<Value<std::vector<int>>>> &
  mIntValues() const {
    return MultiIntValues;
  };

  virtual std::shared_ptr<Value<std::vector<double>>>
  mDoubleValue(size_t i) const {
    return MultiDoubleValues.at(i);
  };

  virtual std::vector<std::shared_ptr<Value<std::vector<double>>>> &
  mDoubleValues() {
    return MultiDoubleValues;
  };

  virtual const std::vector<std::shared_ptr<Value<std::vector<double>>>> &
  mDoubleValues() const {
    return MultiDoubleValues;
  };

  virtual std::shared_ptr<Value<std::vector<std::complex<double>>>>
  mComplexValue(size_t i) const {
    return MultiComplexValues.at(i);
  };

  virtual std::vector<std::shared_ptr<Value<std::vector<std::complex<double>>>>>
      &mComplexValues() {
    return MultiComplexValues;
  };

  virtual const std::vector<
      std::shared_ptr<Value<std::vector<std::complex<double>>>>> &
  mComplexValues() const {
    return MultiComplexValues;
  };

  friend std::ostream &operator<<(std::ostream &out, const ParameterList &b) {
    return out << b.to_str();
  }

  /// A public function returning a string with parameter information
  virtual std::string to_str() const;

protected:
  std::vector<std::shared_ptr<Value<int>>> IntValues;

  std::vector<std::shared_ptr<Value<double>>> DoubleValues;

  std::vector<std::shared_ptr<Value<std::complex<double>>>> ComplexValues;

  std::vector<std::shared_ptr<Value<std::vector<int>>>> MultiIntValues;

  std::vector<std::shared_ptr<Value<std::vector<double>>>> MultiDoubleValues;

  std::vector<std::shared_ptr<Value<std::vector<std::complex<double>>>>>
      MultiComplexValues;

  std::vector<std::shared_ptr<FitParameter>> FitParameters;

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    using namespace boost::serialization;
    // currently only FitParameters can be serialized
    ar &make_nvp("FitParameters", FitParameters);
  }
};

/// Search ParameterList for a FitParameter with \p name. The first match is
/// returned. Be aware that name are not unique. In case no match is found
/// a BadParameter exception is thrown.
inline std::shared_ptr<FitParameter> FindParameter(std::string name,
                                                   const ParameterList &v) {
  auto it =
      std::find_if(v.doubleParameters().begin(), v.doubleParameters().end(),
                   [name](const std::shared_ptr<FitParameter> &s) {
                     return s->name() == name;
                   });
  if (it == v.doubleParameters().end())
    throw BadParameter("FindParameter() | Parameter not in list!");
  return *it;
}

/// Search list for a FitParameter with \p name. The first match is
/// returned. Be aware that name are not unique. In case no match is found
/// a BadParameter exception is thrown.
inline std::shared_ptr<FitParameter>
FindParameter(std::string name, std::vector<std::shared_ptr<FitParameter>> &v) {
  auto it = std::find_if(v.begin(), v.end(),
                         [name](const std::shared_ptr<FitParameter> &s) {
                           return s->name() == name;
                         });
  if (it == v.end())
    throw BadParameter("FindParameter() | Parameter not in list!");
  return *it;
}

inline std::shared_ptr<Value<std::vector<double>>>
findMDoubleValue(const std::string &name, const ParameterList &list) {
  auto it = std::find_if(
      list.mDoubleValues().begin(), list.mDoubleValues().end(),
      [name](const std::shared_ptr<Value<std::vector<double>>> &s) {
        return s->name() == name;
      });
  if (it == list.mDoubleValues().end())
    throw BadParameter("FindParameter() | Parameter not in list!");
  return *it;
}

} // namespace FunctionTree
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

} // namespace serialization
} // namespace boost
#endif // END serialization work-a-round

#endif
