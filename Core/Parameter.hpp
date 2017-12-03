// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Implementations of Parameter for various data types.
///

#ifndef _PARAMETER_HPP_
#define _PARAMETER_HPP_

#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <stdexcept>
#include <cmath>

#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/optional.hpp>

#include "Core/AbsParameter.hpp"
#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"

enum ErrorType { SYM = 1, ASYM = 2, LHSCAN = 3, NOTDEF = 0 };

namespace ComPWA {

///
/// \class MultiComplex
/// Implementation of Parameter for std::vector<std::complex>> => MultiComplex.
///
class MultiComplex : public Parameter {

public:
  /// Standard constructor with just a value provided. Creates parameter
  /// with given value but without bounds or an error.
  MultiComplex(std::string inName,
               const std::vector<std::complex<double>> &values)
      : Parameter(inName, ParType::MCOMPLEX), val_(values) {}

  /// Simple copy constructor using the = operator. As this operator is not
  /// overloaded in this class, c++ will copy every member variable. As this
  /// is a container class, this should be fine.
  MultiComplex(const MultiComplex &in)
      : Parameter(in.Name, ParType::MCOMPLEX) {
    *this = in;
  }

  virtual double numValues() const { return val_.size(); }

  /// Constant reference on values
  virtual const std::vector<std::complex<double>> &values() const {
    return val_;
  }

  virtual std::complex<double> value(unsigned int i) const;

  virtual void setValue(const std::complex<double> inVal, unsigned int i = 0);

  virtual void setValues(std::vector<std::complex<double>> v) { val_ = v; }

  std::vector<std::complex<double>>::iterator first() { return val_.begin(); }

  std::vector<std::complex<double>>::iterator last() { return val_.end(); }

  bool operator==(const MultiComplex otherPar) const;

  virtual std::string to_str() const;

  virtual std::string val_to_str() const;

protected:
  virtual std::string className() const { return "complex collection"; }

  /// Parameter values
  std::vector<std::complex<double>> val_;
};

class MultiDouble : public Parameter {

public:
  /// Standard constructor with a value
  MultiDouble(std::string inName, const std::vector<double> &values)
      : Parameter(inName, ParType::MDOUBLE), val_(values) {}

  /// Copy constructor using = operator.
  /// Simple copy constructor using the = operator. As this operator is not
  /// overloaded in this class, c++ will copy every member variable. As this
  /// is a container class, this should be fine.
  /// \p in input PWAParameter which variables will be copied
  MultiDouble(const MultiDouble &in) : Parameter(in.Name, ParType::MDOUBLE) {
    *this = in;
  }

  virtual double numValues() const { return val_.size(); }

  virtual const std::vector<double> &values() const { return val_; }

  virtual double value(unsigned int i = 0) const;

  virtual void setValue(const double inVal, unsigned int i);

  virtual void setValues(std::vector<double> v) { val_ = v; }

  std::vector<double>::iterator first() { return val_.begin(); }

  std::vector<double>::iterator last() { return val_.end(); }

  bool operator==(const MultiDouble otherPar) const;

  virtual std::string to_str() const;

  virtual std::string val_to_str() const;

protected:
  virtual std::string className() const { return "double collection"; }

  /// Parameter values
  std::vector<double> val_;
};

class MultiUnsignedInteger : public Parameter {

public:
  /// Constructor with just a value provided. Creates parameter
  /// with given value but without bounds or an error.
  MultiUnsignedInteger(std::string inName,
                       const std::vector<unsigned int> &values)
      : Parameter(inName, ParType::MUNSIGNEDINTEGER), val_(values) {}

  /// Simple copy constructor using the = operator. As this operator is not
  /// overloaded in this class, c++ will copy every member variable. As this
  /// is a container class, this should be fine.
  /// param in input PWAParameter which variables will be copied
  MultiUnsignedInteger(const MultiUnsignedInteger &in)
      : Parameter(in.Name, ParType::MUNSIGNEDINTEGER) {
    *this = in;
  }

  virtual double numValues() const { return val_.size(); }

  virtual const std::vector<unsigned int> &values() const { return val_; }

  virtual unsigned int value(unsigned int i = 0) const;

  virtual void setValue(const unsigned int inVal, unsigned int i = 0);

  virtual void setValues(std::vector<unsigned int> v) { val_ = v; };

  bool operator==(const MultiUnsignedInteger otherPar) const;

  virtual std::string to_str() const;

  virtual std::string val_to_str() const;

protected:
  virtual std::string className() const { return "unsigned int collection"; }

  /// Parameter values
  std::vector<unsigned int> val_;
};

class ComplexParameter : public Parameter {

public:
  /// Standard constructor with just a name provided. Creates parameter
  /// with value 0 but without bounds or an error.
  ComplexParameter(std::string inName);

  /// Standard constructor with just value and name provided. Creates parameter
  /// with given value but without bounds or an error.
  ComplexParameter(std::string inName, const std::complex<double> value);

  /// Simple copy constructor using the = operator. As this operator is not
  /// overloaded in this class, c++ will copy every member variable. As this
  /// is a container class, this should be fine.
  ComplexParameter(const ComplexParameter &in);

  /// Parameter value
  virtual std::complex<double> value() const { return val_; }

  /// Set parameter value
  virtual void setValue(const std::complex<double> inVal);

  bool operator==(const ComplexParameter &otherPar) const;

  virtual std::string to_str() const;

  virtual std::string val_to_str() const;

protected:
  virtual std::string className() const { return "complex double"; }

  /// Parameter value
  std::complex<double> val_;
};

//============================================================================
//================================= DOUBLE ===================================
//============================================================================

class DoubleParameter : public Parameter {

public:
  /// Standard constructor with no information provided. Creates parameter
  /// with value 0 but without bounds or an error.
  /// \param inName internal string identifier of this parameter
  DoubleParameter(std::string inName = "");

  /// Construct parameter from a property tree. The expected tree layout is
  /// described in load().
  DoubleParameter(const boost::property_tree::ptree pt);

  /// Standard constructor with just a value provided. Creates parameter
  /// with given value but without bounds or an error.
  DoubleParameter(std::string inName, const double value);

  /// Standard constructor with value and error provided. Creates parameter
  /// with given value and error but without bounds.
  DoubleParameter(std::string inName, const double value, const double error);

  /// Standard constructor with value and bounds provided. Creates parameter
  /// with given value and bounds but without error. If a check for valid
  /// bounds fails, just the value is used.
  DoubleParameter(std::string inName, const double value, const double min,
                  const double max);

  /// Standard constructor with value, bounds and error provided. Creates
  /// parameter with the given information. If a check for valid bounds
  /// fails, just value and error are used.
  DoubleParameter(std::string inName, const double value, const double min,
                  const double max, const double error);

  DoubleParameter(const DoubleParameter &in);

  operator double() const { return Value; };

  virtual bool hasBounds() const { return HasBounds; }

  virtual bool isFixed() const { return IsFixed; }

  virtual void fixParameter(const bool fixed) { IsFixed = fixed; }

  /// Update member variables from other DoubleParameter.
  /// Do to the Observer pattern we can't use a copy constructor.
  /// Therefore we use this workaround. The function ignores if parameter
  /// is fixed!
  virtual void updateParameter(std::shared_ptr<DoubleParameter> newPar);

  //====== PARAMETER VALUE ========
  /// Getter for value of parameter
  virtual double value() const { return Value; }

  /// Setter for value of parameter
  virtual void setValue(const double inVal);

  /// Bounds of parameter
  virtual std::pair<double, double> bounds() const;

  /// Bounds of parameter
  virtual void setBounds(const double min, const double max);

  /// Bounds of parameter
  virtual void setBounds(const std::pair<double, double> r);

  //====== PARAMETER ERROR ========
  /// Is an error set?
  virtual bool hasError() const;

  virtual ErrorType errorType() const { return ErrType; }

  /// Parameter error.
  virtual std::pair<double, double> error() const;

  /// Set parameter error and assume that this parameter has asymmetri errors.
  virtual void setError(double errLow, double errHigh);

  /// Set parameter error and assume that this parameter has asymmetric errors.
  virtual void setError(std::pair<double, double> err);

  /// Setter parameter error and assume that this parameter has symmetric
  /// errors.
  virtual void setError(double err);

  bool operator==(const DoubleParameter otherPar) const;

  /// Load parameters from a ptree. This approach is more or
  /// less equivalent to the serialization of a parameter but provides a better
  /// readable format.
  void load(const boost::property_tree::ptree pt);

  /// Save parameter to a ptree. This approach is more or
  /// less equivalent to the serialization of a parameter but provides a better
  /// readable format.
  boost::property_tree::ptree save() const;

  /// String with detailed information about the parameter. Used in
  /// operator<<().
  virtual std::string to_str() const;

  /// String with detailed information about the parameter. Used in
  /// operator<<().
  virtual std::string val_to_str() const;

protected:
  virtual std::string className() const { return "Double"; }

  /// Are valid bounds defined for this parameter?
  bool HasBounds;

  /// Do you want to keep parameter fixed?
  bool IsFixed;

  /// Parameter value
  double Value;

  /// Parameter bounds
  std::pair<double, double> Bounds;

  /// No error / symmetric error / asymmetric error
  ErrorType ErrType;

  /// Lower parameter error
  std::pair<double, double> Error;

  virtual void SetErrorType(ErrorType t) { ErrType = t; }

  /// Check if \p min and \p max are valid bounds
  bool check_bounds(const std::pair<double, double> bounds) const;

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    using namespace boost::serialization;
    //    ar &boost::serialization::make_nvp(
    //        "Parameter", boost::serialization::base_object<Parameter>(*this));
    ar &make_nvp("Name", Name);
    ar &make_nvp("Bounds", Bounds);
    ar &make_nvp("Fix", IsFixed);
    ar &make_nvp("Value", Value);
    ar &make_nvp("Bounds", Bounds);
    try {
      ar &make_nvp("ErrorType", ErrType);
      ar &make_nvp("Error", Error);
    } catch (...) {
      Error = std::pair<double, double>(0, 0);
      ErrType = ErrorType::SYM;
    }
  }
};
BOOST_SERIALIZATION_SHARED_PTR(ComPWA::DoubleParameter)

class IntegerParameter : public Parameter {

public:
  /// Standard constructor with no information provided. Creates parameter
  /// with value 0 but without bounds or an error.
  /// \param inName internal string identifier of this parameter
  IntegerParameter(std::string inName);

  /// Standard constructor with just a value provided. Creates parameter
  /// with given value but without bounds or an error.
  IntegerParameter(std::string inName, const int value);

  /// Standard constructor with value and error provided. Creates parameter
  /// with given value and error but without bounds.
  IntegerParameter(std::string inName, const int value, const int error);

  /// Standard constructor with value and bounds provided. Creates parameter
  /// with given value and bounds but without error. If a check for valid
  /// bounds fails, just the value is used.
  IntegerParameter(std::string inName, const int value, const int min,
                   const int max);

  /// Standard constructor with value, bounds and error provided. Creates
  /// parameter with the given information. If a check for valid bounds
  /// fails, just value and error are used.
  IntegerParameter(std::string inName, const int value, const int min,
                   const int max, const int error);

  /// Simple copy constructor using the = operator. As this operator is not
  /// overloaded in this class, c++ will copy every member variable. As this
  /// is a container class, this should be fine.
  IntegerParameter(const IntegerParameter &in);

  virtual bool hasBounds() const { return HasBounds; }

  virtual bool hasError() const;

  virtual bool isFixed() const { return IsFixed; }

  virtual int value() const { return Value; }

  virtual std::pair<int, int> bounds() const { return Bounds; }

  virtual std::pair<int, int> error() const { return Error; }

  virtual void setValue(const int inVal);

  virtual void setError(const int inErr);
  
  virtual void setError(const std::pair<int,int> inErr);

  virtual void setBounds(const int min, const int max);

  virtual void setBounds(const std::pair<int, int> r);

  /// Set parameter free or fixed
  virtual void fixParameter(const bool fixed) { IsFixed = fixed; }

  operator int() const { return Value; };

  bool operator==(const IntegerParameter otherPar) const;

  virtual std::string to_str() const;

  virtual std::string val_to_str() const;
  
protected:
  virtual std::string className() const { return "integer"; }

  /// Parameter value
  int Value;
  
  /// Do you want to keep parameter fixed?
  bool IsFixed;
  
  /// Parameter minimum bound
  std::pair<int, int> Bounds;
  
  /// Are valid bounds defined for this parameter?
  bool HasBounds;

  /// No error / symmetric error / asymmetric error
  ErrorType ErrType;

  /// Parameter error
  std::pair<int, int> Error;

 /// This function checks if the bounds of the parameter are valid:
  /// Upper bound should be larger then lower bound and the value
  /// should be inside of the bounds.
  bool check_bounds(const int min, const int max) const;

};

class BoolParameter : public Parameter {

public:
  /// Standard constructor with no information provided. Creates parameter
  /// with value 0 but without bounds or an error.
  /// \param inName internal string identifier of this parameter
  BoolParameter(std::string inName);

 /// Standard constructor with just a value provided. Creates parameter
   /// with given value but without bounds or an error.
  BoolParameter(std::string inName, const bool value);

   /// Simple copy constructor using the = operator. As this operator is not
  /// overloaded in this class, c++ will copy every member variable. As this
  /// is a container class, this should be fine.
  BoolParameter(const BoolParameter &in);

  virtual bool value() const { return Value; }

  virtual void setValue(const bool inVal);

  operator bool() const { return Value; };

  bool operator==(const BoolParameter otherPar) const;

  virtual std::string to_str() const;

  virtual std::string val_to_str() const;
  
protected:
  virtual std::string className() const { return "boolean"; }

  /// Parameter value
  int Value;
};

} // namespace ComPWA

#endif
