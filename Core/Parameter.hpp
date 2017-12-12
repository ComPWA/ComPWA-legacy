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

namespace ComPWA {
enum ErrorType { SYM = 1, ASYM = 2, LHSCAN = 3, NOTDEF = 0 };

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
  
  virtual bool isParameter() const { return true; }

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

  /// Average arameter error (in case of asymmetric errors) or simpli parameter
  /// error.
  virtual double avgError() const { return 0.5 * (Error.first + Error.second); }

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

} // namespace ComPWA

#endif
