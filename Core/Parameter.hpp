// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Implementations for internal parameter.
/*! \class DoubleParameter
 * \class IntegerParameter
 * \class BoolParameter
 * @file Parameter.hpp
 * This class implements some internal container of parameters.
 * A parameter consists of a value with optional bounds and error.
 */

#ifndef _PARAMETER_HPP_
#define _PARAMETER_HPP_

#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <stdexcept>
#include <cmath>

#include <boost/property_tree/ptree.hpp>

#include "Core/AbsParameter.hpp"
#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"

enum ErrorType { SYM = 1, ASYM = 2, LHSCAN = 3, NOTDEF = 0 };

namespace ComPWA {

class MultiComplex : public AbsParameter {

public:
  //! Standard constructor without information
  /*!
   * Standard constructor with no information provided. Creates parameter
   * with value 0 but without bounds or an error.
   * \param inName internal string identifier of this parameter
   */
  // MultiComplex(std::string inName):AbsParameter(inName, ParType::MDOUBLE){
  //}
  //! Standard constructor with a value
  /*!
   * Standard constructor with just a value provided. Creates parameter
   * with given value but without bounds or an error.
   * \param inName internal string identifier of this parameter
   * \param values input vector of values of the parameter
   */
  MultiComplex(std::string inName,
               const std::vector<std::complex<double>> &values)
      : AbsParameter(inName, ParType::MCOMPLEX), val_(values) {}

  //! Copy constructor using = operator
  /*!
   * Simple copy constructor using the = operator. As this operator is not
   * overloaded in this class, c++ will copy every member variable. As this
   * is a container class, this should be fine.
   * \param in input PWAParameter which variables will be copied
   */
  MultiComplex(const MultiComplex &in)
      : AbsParameter(in.name_, ParType::MCOMPLEX) {
    *this = in;
    //      error_ = std::shared_ptr<ParError<double>>(new
    //      ParError<double>(*in.error_));
  }

  //! Getter for number of values in this multipar
  virtual double GetNValues() const { return val_.size(); }

  //! Getter for value of parameter
  virtual const std::vector<std::complex<double>> &GetValues() const {
    return val_;
  }

  virtual std::complex<double> GetValue(unsigned int i = 0) const;

  virtual void SetValue(const std::complex<double> inVal, unsigned int i = 0);

  std::vector<std::complex<double>>::iterator Begin() { return val_.begin(); }

  std::vector<std::complex<double>>::iterator End() { return val_.end(); }

  bool operator==(const MultiComplex otherPar) const;

protected:
  virtual std::string TypeName() const { return "complex collection"; }

  /// Parameter values
  std::vector<std::complex<double>> val_;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses all available information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str(), type()
   */
  virtual std::string make_str() const;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses only the value information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, make_str()
   */
  virtual std::string make_val_str() const;
};

class MultiDouble : public AbsParameter {

public:
  //! Standard constructor without information
  /*!
   * Standard constructor with no information provided. Creates parameter
   * with value 0 but without bounds or an error.
   * \param inName internal string identifier of this parameter
   */
  // MultiDouble(std::string inName):AbsParameter(inName, ParType::MDOUBLE){
  //}
  //! Standard constructor with a value
  /*!
   * Standard constructor with just a value provided. Creates parameter
   * with given value but without bounds or an error.
   * \param inName internal string identifier of this parameter
   * \param values input vector of values of the parameter
   */
  MultiDouble(std::string inName, const std::vector<double> &values)
      : AbsParameter(inName, ParType::MDOUBLE), val_(values) {}

  //! Copy constructor using = operator
  /*!
   * Simple copy constructor using the = operator. As this operator is not
   * overloaded in this class, c++ will copy every member variable. As this
   * is a container class, this should be fine.
   * \param in input PWAParameter which variables will be copied
   */
  MultiDouble(const MultiDouble &in)
      : AbsParameter(in.name_, ParType::MDOUBLE) {
    *this = in;
    //		error_ = std::shared_ptr<ParError<double>>(new
    // ParError<double>(*in.error_));
  }

  //! Get number of values in MultiDouble parameter
  virtual double GetNValues() const { return val_.size(); }

  virtual const std::vector<double> &GetValues() const { return val_; }

  virtual double GetValue(unsigned int i = 0) const;

  virtual void SetValue(const double inVal, unsigned int i = 0);

  std::vector<double>::iterator Begin() { return val_.begin(); }

  std::vector<double>::iterator End() { return val_.end(); }

  bool operator==(const MultiDouble otherPar) const;

protected:
  virtual std::string TypeName() const { return "double collection"; }

  /// Parameter values
  std::vector<double> val_;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses all available information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str(), type()
   */
  virtual std::string make_str() const;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses only the value information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, make_str()
   */
  virtual std::string make_val_str() const;
};

class MultiUnsignedInteger : public AbsParameter {

public:
  //! Standard constructor with a value
  /*!
   * Standard constructor with just a value provided. Creates parameter
   * with given value but without bounds or an error.
   * \param inName internal string identifier of this parameter
   * \param values input vector of values of the parameter
   */
  MultiUnsignedInteger(std::string inName,
                       const std::vector<unsigned int> &values)
      : AbsParameter(inName, ParType::MUNSIGNEDINTEGER), val_(values) {}

  //! Copy constructor using = operator
  /*!
   * Simple copy constructor using the = operator. As this operator is not
   * overloaded in this class, c++ will copy every member variable. As this
   * is a container class, this should be fine.
   * \param in input PWAParameter which variables will be copied
   */
  MultiUnsignedInteger(const MultiUnsignedInteger &in)
      : AbsParameter(in.name_, ParType::MUNSIGNEDINTEGER) {
    *this = in;
    //    error_ = std::shared_ptr<ParError<double>>(new
    //    ParError<double>(*in.error_));
  }

  //! Getter for number of values in this multipar
  virtual double GetNValues() const { return val_.size(); }

  virtual const std::vector<unsigned int> &GetValues() const { return val_; }

  virtual unsigned int GetValue(unsigned int i = 0) const;

  virtual void SetValue(const unsigned int inVal, unsigned int i = 0);

  bool operator==(const MultiUnsignedInteger otherPar) const;

protected:
  virtual std::string TypeName() const { return "unsigned int collection"; }

  /// Parameter values
  std::vector<unsigned int> val_;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses all available information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str(), type()
   */
  virtual std::string make_str() const;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses only the value information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, make_str()
   */
  virtual std::string make_val_str() const;

private:
};

class ComplexParameter : public AbsParameter {

public:
  //! Standard constructor without information
  /*!
   * Standard constructor with just a name provided. Creates parameter
   * with value 0 but without bounds or an error.
   * \param inName internal string identifier of this parameter
   */
  ComplexParameter(std::string inName);

  //! Standard constructor with a value
  /*!
   * Standard constructor with just value and name provided. Creates parameter
   * with given value but without bounds or an error.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   */
  ComplexParameter(std::string inName, const std::complex<double> value);

  //! Standard constructor with value and error
  /*!
   * Standard constructor with value and error provided. Creates parameter
   * with given value and error but without bounds.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   * \param error input error of the parameter
   */
  ComplexParameter(std::string inName, const std::complex<double> value,
                   const std::complex<double> error);

  //! Standard constructor with value and bounds
  /*!
   * Standard constructor with value and bounds provided. Creates parameter
   * with given value and bounds but without error. If a check for valid
   * bounds fails, just the value is used.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   * \param min input lower bound
   * \param max input upper bound
   * \sa check_bounds()
   */
  ComplexParameter(std::string inName, const std::complex<double> value,
                   const std::complex<double> min,
                   const std::complex<double> max);

  //! Standard constructor with value, bounds and error
  /*!
   * Standard constructor with value, bounds and error provided. Creates
   * parameter with the given information. If a check for valid bounds
   * fails, just value and error are used.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   * \param min input lower bound
   * \param max input upper bound
   * \param error input error of the parameter
   * \sa check_bounds()
   */
  ComplexParameter(std::string inName, const std::complex<double> value,
                   const std::complex<double> min,
                   const std::complex<double> max,
                   const std::complex<double> error);

  //! Copy constructor using = operator
  /*!
   * Simple copy constructor using the = operator. As this operator is not
   * overloaded in this class, c++ will copy every member variable. As this
   * is a container class, this should be fine.
   * \param in input PWAParameter which variables will be copied
   */
  ComplexParameter(const ComplexParameter &in);

  //! Check if parameter has bounds
  virtual bool HasBounds() const { return bounds_; }

  //! Check if bounds should be used
  virtual bool UseBounds() const {
    if (bounds_)
      return usebounds_;
    return false;
  }

  //! Check if parameter has an error
  virtual bool HasError() const { return hasError_; }

  //! Check if parameter is fixed
  virtual bool IsFixed() const { return fixed_; }

  //! Getter for value of parameter
  virtual std::complex<double> GetValue() const { return val_; }

  //! Getter for lower bound of parameter
  virtual std::complex<double> GetMinValue() const { return min_; }

  //! Getter for upper bound of parameter
  virtual std::complex<double> GetMaxValue() const { return max_; }

  //! Getter for error of parameter
  virtual std::complex<double> GetError() const { return err_; }

  void SetError(const std::complex<double> inErr);

  //! Setter for value of parameter
  virtual void SetValue(const std::complex<double> inVal);

  //! Setter for lower bound
  /*!
   * Setter for lower bound of the parameter. If a check for valid bounds
   * fails, it returns false and nothing changes. This means if the lower
   * bound is invalid the parameter maintains its old bounds if it had some.
   * \param min input lower bound
   * \return bool if successful (re)set lower bound
   * \sa check_bounds()
   */
  virtual bool SetMinValue(const std::complex<double> min);

  //! Setter for upper bound
  /*!
   * Setter for upper bound of the parameter. If a check for valid bounds
   * fails, it returns false and nothing changes. This means if the upper
   * bound is invalid the parameter maintains its old bounds if it had some.
   * \param max input upper bound
   * \return bool if successful (re)set upper bound
   * \sa check_bounds()
   */
  virtual bool SetMaxValue(const std::complex<double> max);

  virtual bool SetMinMax(const std::complex<double> inMin,
                         const std::complex<double> inMax);

  //! Set if bounds should be used
  virtual void UseBounds(const bool use) { usebounds_ = use; }

  //! Call to fix parameter
  virtual void SetParameterFixed() { fixed_ = true; }

  //! Call to free parameter
  virtual void SetParameterFree() { fixed_ = false; }

  //! Set parameter free or fixed
  virtual void FixParameter(const bool fixed) { fixed_ = fixed; }

  bool operator==(const ComplexParameter &otherPar) const;

protected:
  virtual std::string TypeName() const { return "complex double"; }

  bool bounds_; /*!< Are valid bounds defined for this parameter? */

  bool hasError_; /*!< Is an error defined for this parameter? */

  bool usebounds_; /*!< Do you want to restrict your parameter? */

  bool fixed_; /*!< Do you want to keep parameter fixed? */

  /// Parameter value
  std::complex<double> val_;

  /// Parameter error
  std::complex<double> err_;

  /// Parameter minimum bound
  std::complex<double> min_;

  /// Parameter maximum bound
  std::complex<double> max_;

  //! A protected function to check if bounds are valid
  /*!
   * This function checks if the bounds of the parameter are valid:
   * Upper bound should be larger then lower bound and the value
   * should be inside of the bounds.
   * \param max upper bound to check
   * \param min lower bound to check
   * \return bool if bounds are valid
   * \sa Parameter(const double value, const double min, const double max)
   * \sa Parameter(const double value, const double min, const double max, const
   * double error)
   * \sa SetMinMax(), SetMinValue(), SetMaxValue()
   */
  bool check_bounds(const std::complex<double> min,
                    const std::complex<double> max);

  //! A protected function which returns an output string for printing
  /*!
   * This function uses all available information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str(), type()
   */
  virtual std::string make_str() const;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses only the value information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, make_str()
   */
  virtual std::string make_val_str() const;
};

//============================================================================
//================================= DOUBLE ===================================
//============================================================================

class DoubleParameter : public AbsParameter {

public:
  //! Standard constructor without information
  /*!
   * Standard constructor with no information provided. Creates parameter
   * with value 0 but without bounds or an error.
   * \param inName internal string identifier of this parameter
   */
  DoubleParameter(std::string inName = "");

  //! Standard constructor with a value
  /*!
   * Standard constructor with just a value provided. Creates parameter
   * with given value but without bounds or an error.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   */
  DoubleParameter(std::string inName, const double value);

  //! Standard constructor with value and error
  /*!
   * Standard constructor with value and error provided. Creates parameter
   * with given value and error but without bounds.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   * \param error input error of the parameter
   */
  DoubleParameter(std::string inName, const double value, const double error);

  //! Standard constructor with value and bounds
  /*!
   * Standard constructor with value and bounds provided. Creates parameter
   * with given value and bounds but without error. If a check for valid
   * bounds fails, just the value is used.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   * \param min input lower bound
   * \param max input upper bound
   * \sa check_bounds()
   */
  DoubleParameter(std::string inName, const double value, const double min,
                  const double max);

  //! Standard constructor with value, bounds and error
  /*!
   * Standard constructor with value, bounds and error provided. Creates
   * parameter with the given information. If a check for valid bounds
   * fails, just value and error are used.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   * \param min input lower bound
   * \param max input upper bound
   * \param error input error of the parameter
   * \sa check_bounds()
   */
  DoubleParameter(std::string inName, const double value, const double min,
                  const double max, const double error);

  DoubleParameter(const DoubleParameter &in);

  //! Operator for conversion to double
  operator double() const { return val_; };

  //! Check if parameter has bounds
  virtual bool HasBounds() const { return bounds_; }

  //! Check if parameter is fixed
  virtual bool IsFixed() const { return fixed_; }

  //! Call to fix parameter
  virtual void SetParameterFixed() { fixed_ = true; }

  //! Call to free parameter
  virtual void SetParameterFree() { fixed_ = false; }

  //! Set parameter free or fixed
  virtual void FixParameter(const bool fixed) { fixed_ = fixed; }

  /*! Update member variables from other DoubleParameter
   * Do to the Observer pattern we can't use a copy constructor.
   * Therefore we use this workaround. The function ignores if parameter
   * is fixed!
   */
  virtual void UpdateParameter(std::shared_ptr<DoubleParameter> newPar);

  //====== PARAMETER VALUE ========
  //! Getter for value of parameter
  virtual double GetValue() const { return val_; }

  //! Getter for value of parameter
  virtual double GetRoundedValue() const { return val_; }

  //! Getter for lower bound of parameter
  virtual double GetMinValue() const { return min_; }

  //! Getter for upper bound of parameter
  virtual double GetMaxValue() const { return max_; }

  //! Setter for value of parameter
  virtual void SetValue(const double inVal);

  //! Setter for bounds of parameter
  virtual void SetRange(const double min, const double max);

  //! Setter for bounds of parameter
  virtual void SetMinMax(const double min, const double max);

  /*! Setter for lower bound
   * Setter for lower bound of the parameter. If a check for valid bounds
   * fails, it returns false and nothing changes. This means if the lower
   * bound is invalid the parameter maintains its old bounds if it had some.
   * \param min input lower bound
   * \return bool if successful (re)set lower bound
   * \sa check_bounds()
   */
  virtual void SetMinValue(const double min);

  /*! Setter for upper bound
   * Setter for upper bound of the parameter. If a check for valid bounds
   * fails, it returns false and nothing changes. This means if the upper
   * bound is invalid the parameter maintains its old bounds if it had some.
   * \param max input upper bound
   * \return bool if successful (re)set upper bound
   * \sa check_bounds()
   */
  virtual void SetMaxValue(const double max);

  //====== PARAMETER ERROR ========
  //! Check if parameter has an error
  virtual bool HasError() const;

  //! Getter for type of parameter error
  virtual ErrorType GetErrorType() const { return errorType; }

  //! Getter for parameter error. In case of asymmetric errors the average error
  //! is returned.
  virtual double GetError() const;

  //! Get rounded parameter error. In case of asymmetric errors the average
  //! error is returned.
  virtual double GetRoundedError() const { return GetError(); }

  //! Getter for upper error of parameter
  virtual double GetErrorHigh() const;

  //! Getter for lower error of parameter
  virtual double GetErrorLow() const;

  //! Setter for low/high error of parameter
  virtual void SetError(double errLow, double errHigh);

  //! Setter for error of parameter
  virtual void SetError(double err);

  bool operator==(const DoubleParameter otherPar) const;

protected:
  virtual std::string TypeName() const { return "double"; }

  bool bounds_; /*!< Are valid bounds defined for this parameter? */

  bool fixed_; /*!< Do you want to keep parameter fixed? */

  /// Parameter value
  double val_;

  /// Parameter lower bound
  double min_;

  /// Parameter upper bound
  double max_;

  /// error type
  ErrorType errorType;

  /// Lower parameter error
  double errorLow;

  /// Upper parameter error
  double errorHigh;

  virtual void SetErrorHigh(double errHigh) { errorHigh = errHigh; }

  virtual void SetErrorLow(double errLow) { errorLow = std::fabs(errLow); }

  virtual void SetErrorType(ErrorType t) { errorType = t; }

  //! A protected function to check if bounds are valid
  /*!
   * This function checks if the bounds of the parameter are valid:
   * Upper bound should be larger then lower bound and the value
   * should be inside of the bounds.
   * \param max upper bound to check
   * \param min lower bound to check
   * \return bool if bounds are valid
   * \sa Parameter(const double value, const double min, const double max)
   * \sa Parameter(const double value, const double min, const double max, const
   * double error)
   * \sa SetMinMax(), SetMinValue(), SetMaxValue()
   */
  bool check_bounds(const double min, const double max);

  //! A protected function which returns an output string for printing
  /*!
   * This function uses all available information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str(), type()
   */
  virtual std::string make_str() const;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses only the value information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, make_str()
   */
  virtual std::string make_val_str() const;

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    using namespace boost::serialization;
    ar &boost::serialization::make_nvp(
        "AbsParameter", boost::serialization::base_object<AbsParameter>(*this));
    ar &make_nvp("bounds", bounds_);
    ar &make_nvp("isFixed", fixed_);
    ar &make_nvp("value", val_);
    ar &make_nvp("min_value", min_);
    ar &make_nvp("max_value", max_);
    try {
      ar &make_nvp("errorType", errorType);
      ar &make_nvp("errorLow", errorLow);
      ar &make_nvp("errorHigh", errorHigh);
    } catch (...) {
      errorLow = 0;
      errorHigh = 0;
      errorType = ErrorType::SYM;
    }
  }
};
BOOST_SERIALIZATION_SHARED_PTR(ComPWA::DoubleParameter)

/**
 Create a DoubleParameter object from a ptree. This approach is more or
 less equivalent to the serialization of a parameter but provides a better
 readable format.

 @param pt Input property tree
 @return Parameter
 */
inline DoubleParameter
DoubleParameterFactory(const boost::property_tree::ptree pt) {
  DoubleParameter obj;

  // Require that name and value are provided
  obj.SetName(pt.get<std::string>("<xmlattr>.Name"));
  obj.SetValue(pt.get<double>("Value"));

  // Optional settings
  if (pt.get_optional<double>("Error")) {
    obj.SetError(pt.get<double>("Error"));
  }
  if (pt.get_optional<double>("ErrorLow")) {
    if (pt.get_optional<double>("ErrorHigh"))
      obj.SetError(pt.get<double>("ErrorLow"), pt.get<double>("ErrorHigh"));
    else
      throw std::runtime_error("DoubleParameterFactory() | Parameter asymmetic "
                               "error not properly set!");
  } else if (pt.get_optional<double>("ErrorHigh")) {
    throw std::runtime_error("DoubleParameterFactory() | Parameter asymmetic "
                             "error not properly set!");
  } else { /* Do not set a asymmetric errors */
  }

  if (pt.get_optional<bool>("Fix"))
    obj.FixParameter(pt.get<bool>("Fix"));

  if (pt.get_optional<double>("Min")) {
    if (pt.get_optional<double>("Max")) {
      obj.SetRange(pt.get<double>("Min"), pt.get<double>("Max"));
    } else {
      throw std::runtime_error(
          "DoubleParameterFactory() | Parameter range not properly set!");
    }
  } else if (pt.get_optional<double>("Max")) {
    throw std::runtime_error(
        "DoubleParameterFactory() | Parameter range not properly set!");
  } else {
    // Do not set a parameter range and fix the parameter
    obj.FixParameter(true);
  }

  return obj;
}

/// Save a DoubleParameter object from a ptree. This approach is more or
/// less equivalent to the serialization of a parameter but provides a better
/// readable format.
inline boost::property_tree::ptree
DoubleParameterSave(const DoubleParameter par) {
  boost::property_tree::ptree pt;

  // Require that name and value are provided
  pt.put("<xmlattr>.Name", par.GetName());
  pt.put("Value", par.GetValue());
  pt.put("Fix", par.IsFixed());
  if (par.HasBounds()) {
    pt.put("Min", par.GetMinValue());
    pt.put("Max", par.GetMaxValue());
  }
  if (par.HasError()) {
    if (par.GetErrorLow() == par.GetErrorHigh()) {
      pt.put("Error", par.GetError());
    } else {
      pt.put("ErrorLow", par.GetErrorLow());
      pt.put("ErrorHigh", par.GetErrorHigh());
    }
  }

  return pt;
}

inline boost::property_tree::ptree
DoubleParameterSave(std::shared_ptr<DoubleParameter> par) {
  return DoubleParameterSave(*par.get());
}

class IntegerParameter : public AbsParameter {

public:
  //! Standard constructor without information
  /*!
   * Standard constructor with no information provided. Creates parameter
   * with value 0 but without bounds or an error.
   * \param inName internal string identifier of this parameter
   */
  IntegerParameter(std::string inName);

  //! Standard constructor with a value
  /*!
   * Standard constructor with just a value provided. Creates parameter
   * with given value but without bounds or an error.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   */
  IntegerParameter(std::string inName, const int value);

  //! Standard constructor with value and error
  /*!
   * Standard constructor with value and error provided. Creates parameter
   * with given value and error but without bounds.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   * \param error input error of the parameter
   */
  IntegerParameter(std::string inName, const int value, const int error);

  //! Standard constructor with value and bounds
  /*!
   * Standard constructor with value and bounds provided. Creates parameter
   * with given value and bounds but without error. If a check for valid
   * bounds fails, just the value is used.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   * \param min input lower bound
   * \param max input upper bound
   * \sa check_bounds()
   */
  IntegerParameter(std::string inName, const int value, const int min,
                   const int max);

  //! Standard constructor with value, bounds and error
  /*!
   * Standard constructor with value, bounds and error provided. Creates
   * parameter with the given information. If a check for valid bounds
   * fails, just value and error are used.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   * \param min input lower bound
   * \param max input upper bound
   * \param error input error of the parameter
   * \sa check_bounds()
   */
  IntegerParameter(std::string inName, const int value, const int min,
                   const int max, const int error);

  //! Copy constructor using = operator
  /*!
   * Simple copy constructor using the = operator. As this operator is not
   * overloaded in this class, c++ will copy every member variable. As this
   * is a container class, this should be fine.
   * \param in input PWAParameter which variables will be copied
   */
  IntegerParameter(const IntegerParameter &in);

  //! Check if parameter has bounds
  virtual bool HasBounds() const { return bounds_; }

  //! Check if bounds should be used
  virtual bool UseBounds() const {
    if (bounds_)
      return usebounds_;
    return false;
  }

  //! Check if parameter has an error
  virtual bool HasError() const { return hasError_; }

  //! Check if parameter is fixed
  virtual bool IsFixed() const { return fixed_; }

  //! Getter for value of parameter
  virtual int GetValue() const { return val_; }

  //! Getter for lower bound of parameter
  virtual int GetMinValue() const { return min_; }

  //! Getter for upper bound of parameter
  virtual int GetMaxValue() const { return max_; }

  //! Getter for error of parameter
  virtual int GetError() const { return err_; }

  //! Setter for value of parameter
  virtual void SetValue(const int inVal);

  //! Setter for error of parameter
  virtual void SetError(const int inErr);

  //! Setter for bounds of parameter
  virtual bool SetMinMax(const int inMin, const int inMax);

  //! Setter for lower bound
  /*!
   * Setter for lower bound of the parameter. If a check for valid bounds
   * fails, it returns false and nothing changes. This means if the lower
   * bound is invalid the parameter maintains its old bounds if it had some.
   * \param min input lower bound
   * \return bool if successful (re)set lower bound
   * \sa check_bounds()
   */
  virtual bool SetMinValue(const int min);

  //! Setter for upper bound
  /*!
   * Setter for upper bound of the parameter. If a check for valid bounds
   * fails, it returns false and nothing changes. This means if the upper
   * bound is invalid the parameter maintains its old bounds if it had some.
   * \param max input upper bound
   * \return bool if successful (re)set upper bound
   * \sa check_bounds()
   */
  virtual bool SetMaxValue(const int max);

  //! Set if bounds should be used
  virtual void UseBounds(const bool use) { usebounds_ = use; }

  //! Call to fix parameter
  virtual void SetParameterFixed() { fixed_ = true; }

  //! Call to free parameter
  virtual void SetParameterFree() { fixed_ = false; }

  //! Set parameter free or fixed
  virtual void FixParameter(const bool fixed) { fixed_ = fixed; }

  //! A public function returning a string naming its type
  /*!
   * This function is used to get the type of the implementation of this
   * general parameter interface.
   * \sa operator<<, to_str(), make_str()
   */
  operator int() const { return val_; };

  bool operator==(const IntegerParameter otherPar) const;

protected:
  virtual std::string TypeName() const { return "integer"; }

  bool bounds_; /*!< Are valid bounds defined for this parameter? */

  bool hasError_; /*!< Is an error defined for this parameter? */

  bool usebounds_; /*!< Do you want to restrict your parameter? */

  bool fixed_; /*!< Do you want to keep parameter fixed? */

  /// Parameter value
  int val_;

  /// Parameter minimum bound
  int min_;

  /// Parameter maximum bound
  int max_;

  /// Parameter error
  int err_;

  //! A protected function to check if bounds are valid
  /*!
   * This function checks if the bounds of the parameter are valid:
   * Upper bound should be larger then lower bound and the value
   * should be inside of the bounds.
   * \param max upper bound to check
   * \param min lower bound to check
   * \return bool if bounds are valid
   * \sa Parameter(const double value, const double min, const double max)
   * \sa Parameter(const double value, const double min, const double max, const
   * double error)
   * \sa SetMinMax(), SetMinValue(), SetMaxValue()
   */
  bool check_bounds(const int min, const int max);

  //! A protected function which returns an output string for printing
  /*!
   * This function uses all available information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str(), type()
   */
  virtual std::string make_str() const;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses only the value information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, make_str()
   */
  virtual std::string make_val_str() const;
};

class BoolParameter : public AbsParameter {

public:
  //! Standard constructor without information
  /*!
   * Standard constructor with no information provided. Creates parameter
   * with value 0 but without bounds or an error.
   * \param inName internal string identifier of this parameter
   */
  BoolParameter(std::string inName);

  //! Standard constructor with a value
  /*!
   * Standard constructor with just a value provided. Creates parameter
   * with given value but without bounds or an error.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   */
  BoolParameter(std::string inName, const bool value);

  //! Standard constructor with value and error
  /*!
   * Standard constructor with value and error provided. Creates parameter
   * with given value and error but without bounds.
   * \param inName internal string identifier of this parameter
   * \param value input value of the parameter
   * \param error input error of the parameter
   */
  BoolParameter(std::string inName, const bool value, const bool error);

  //! Copy constructor using = operator
  /*!
   * Simple copy constructor using the = operator. As this operator is not
   * overloaded in this class, c++ will copy every member variable. As this
   * is a container class, this should be fine.
   * \param in input PWAParameter which variables will be copied
   */
  BoolParameter(const BoolParameter &in);

  //! Check if parameter has an error
  virtual bool HasError() const { return hasError_; }

  //! Check if parameter is fixed
  virtual bool IsFixed() const { return fixed_; }

  //! Getter for value of parameter
  virtual bool GetValue() const { return val_; }

  //! Getter for error of parameter
  virtual bool GetError() const { return err_; }

  //! Setter for value of parameter
  virtual void SetValue(const bool inVal);

  //! Setter for error of parameter
  virtual void SetError(const bool inErr);

  //! Call to fix parameter
  virtual void SetParameterFixed() { fixed_ = true; }

  //! Call to free parameter
  virtual void SetParameterFree() { fixed_ = false; }

  //! Set parameter free or fixed
  virtual void FixParameter(const bool fixed) { fixed_ = fixed; }

  bool operator==(const BoolParameter otherPar) const;

protected:
  virtual std::string TypeName() const { return "boolean"; }

  bool hasError_; /*!< Is an error defined for this parameter? */

  bool usebounds_; /*!< Do you want to restrict your parameter? */

  bool fixed_; /*!< Do you want to keep parameter fixed? */

  /// Parameter value
  int val_;

  /// Parameter error
  int err_;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses all available information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str(), type()
   */
  virtual std::string make_str() const;

  //! A protected function which returns an output string for printing
  /*!
   * This function uses only the value information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, make_str()
   */
  virtual std::string make_val_str() const;
};

} /* namespace ComPWA */

#endif
