
// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Parameter.hpp"

using namespace ComPWA;

//=============================== MULTICOMPEX ================================
//! Getter for value of parameter
std::complex<double> MultiComplex::value(unsigned int i) const {
  if (i >= val_.size())
    return 0;
  return val_.at(i);
}

//! Setter for value of parameter
void MultiComplex::setValue(const std::complex<double> inVal, unsigned int i) {
  if (i >= val_.size())
    return;
  if (val_.at(i) == inVal)
    return;
  val_.at(i) = inVal;
  Notify();
}

bool MultiComplex::operator==(const MultiComplex otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->name() != otherPar.name())
    return false;
  if (this->values() != otherPar.values())
    return false;

  return true;
}

std::string MultiComplex::to_str() const {
  std::stringstream oss;
  oss << Name;
  size_t max = val_.size();
  if (max > 5)
    max = 5; // display only 5 variables
  oss << "\t Val = ";
  for (unsigned int i = 0; i < max - 1; i++)
    oss << val_[i] << ", ";
  oss << val_[max - 1];
  if (max < val_.size())
    oss << " ... ";
  oss << "\t Type = " << className();
  return oss.str();
}

std::string MultiComplex::val_to_str() const {
  std::stringstream ovs;
  size_t max = val_.size();
  if (max > 0) {
    if (max > 3)
      max = 3; // display only 10 variables
    for (unsigned int i = 0; i < max - 1; i++)
      ovs << val_[i] << ", ";
    ovs << val_[max - 1];
    if (max < val_.size())
      ovs << " ... ";
  }
  return ovs.str();
}
//=============================== MULTIDOUBLE ================================

double MultiDouble::value(unsigned int i) const {
  if (i >= val_.size())
    return 0;
  return val_[i];
}

void MultiDouble::setValue(const double inVal, unsigned int i) {
  if (i >= val_.size())
    return;
  if (val_[i] == inVal)
    return;
  val_[i] = inVal;
  Notify();
}

bool MultiDouble::operator==(const MultiDouble otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->name() != otherPar.name())
    return false;
  if (this->value() != otherPar.value())
    return false;

  return true;
}

std::string MultiDouble::to_str() const {
  std::stringstream oss;
  oss << Name;
  size_t max = val_.size();
  if (max > 5)
    max = 5; // display only 10 variables
  oss << "\t Val = ";
  for (unsigned int i = 0; i < max - 1; i++)
    oss << val_[i] << ", ";
  oss << val_[max - 1];
  if (max < val_.size())
    oss << " ... ";
  oss << "\t Type = " << className();
  return oss.str();
}

std::string MultiDouble::val_to_str() const {
  std::stringstream ovs;
  size_t max = val_.size();
  if (max > 0) {
    if (max > 5)
      max = 5; // display only 5 variables
    for (unsigned int i = 0; i < max - 1; i++)
      ovs << val_[i] << ", ";
    ovs << val_[max - 1];
    if (max < val_.size())
      ovs << " ... ";
  }
  return ovs.str();
}
//============================= MULTIUNSIGNEDINT =============================

//! Getter for value of parameter
unsigned int MultiUnsignedInteger::value(unsigned int i) const {
  if (i >= val_.size())
    return 0;
  return val_[i];
}

//! Setter for value of parameter
void MultiUnsignedInteger::setValue(const unsigned int inVal, unsigned int i) {
  if (i >= val_.size())
    return;
  if (val_[i] == inVal)
    return;
  val_[i] = inVal;
  Notify();
}

bool MultiUnsignedInteger::
operator==(const MultiUnsignedInteger otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->name() != otherPar.name())
    return false;
  if (this->value() != otherPar.value())
    return false;

  return true;
}
std::string MultiUnsignedInteger::to_str() const {
  std::stringstream oss;
  oss << Name;
  size_t max = val_.size();
  if (max > 5)
    max = 5; // display only 10 variables
  oss << "\t Val = ";
  for (unsigned int i = 0; i < max - 1; i++)
    oss << val_[i] << ", ";
  oss << val_[max - 1];
  if (max < val_.size())
    oss << " ... ";
  oss << "\t Type = " << className();
  return oss.str();
}

std::string MultiUnsignedInteger::val_to_str() const {
  std::stringstream ovs;
  size_t max = val_.size();
  if (max > 5)
    max = 5; // display only 5 variables
  for (unsigned int i = 0; i < max - 1; i++)
    ovs << val_[i] << ", ";
  ovs << val_[max - 1];
  if (max < val_.size())
    ovs << " ... ";
  return ovs.str();
}
//================================= COMPLEX ==================================
ComplexParameter::ComplexParameter(std::string inName)
    : Parameter(inName, ParType::COMPLEX), val_(0., 0.) {}

ComplexParameter::ComplexParameter(std::string inName,
                                   const std::complex<double> value)
    : Parameter(inName, ParType::COMPLEX), val_(value) {}

ComplexParameter::ComplexParameter(const ComplexParameter &in)
    : Parameter(in.Name, ParType::COMPLEX) {
  *this = in;
}

void ComplexParameter::setValue(const std::complex<double> inVal) {
  if (val_ == inVal)
    return;
  val_ = inVal;
  Notify();
}

bool ComplexParameter::operator==(const ComplexParameter &otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->name() != otherPar.name())
    return false;
  if (this->value() != otherPar.value())
    return false;
  return true;
}

std::string ComplexParameter::to_str() const {
  std::stringstream oss;
  oss << Name;
  oss << "\t Val = " << val_;
  oss << "\t Type = " << className();
  return oss.str();
}

std::string ComplexParameter::val_to_str() const {
  std::stringstream ovs;
  ovs << val_;
  return ovs.str();
}
//================================= DOUBLE ===================================

DoubleParameter::DoubleParameter(std::string inName)
    : Parameter(inName, ParType::DOUBLE), HasBounds(false), IsFixed(0),
      Value(0), Bounds(std::pair<double, double>(0, 0)),
      ErrType(ErrorType::NOTDEF), Error(std::pair<double, double>(0, 0)) {}

DoubleParameter::DoubleParameter(const boost::property_tree::ptree pt)
    : Parameter("", ParType::DOUBLE) {
  load(pt);
}

DoubleParameter::DoubleParameter(std::string inName, const double value)
    : Parameter(inName, ParType::DOUBLE), HasBounds(false), IsFixed(0),
      Value(value), Bounds(std::pair<double, double>(0, 0)),
      ErrType(ErrorType::NOTDEF), Error(std::pair<double, double>(0, 0)) {}

DoubleParameter::DoubleParameter(std::string inName, const double value,
                                 const double error)
    : Parameter(inName, ParType::DOUBLE), HasBounds(false), IsFixed(0),
      Value(value), Bounds(std::pair<double, double>(0, 0)),
      ErrType(ErrorType::SYM), Error(std::pair<double, double>(error, error)) {}

DoubleParameter::DoubleParameter(std::string inName, const double value,
                                 const double min, const double max)
    : Parameter(inName, ParType::DOUBLE), HasBounds(false), IsFixed(0),
      Value(value), Bounds(std::pair<double, double>(0, 0)),
      ErrType(ErrorType::NOTDEF), Error(std::pair<double, double>(0, 0)) {
  setBounds(min, max);
}

DoubleParameter::DoubleParameter(std::string inName, const double value,
                                 const double min, const double max,
                                 const double error)
    : Parameter(inName, ParType::DOUBLE), HasBounds(false), IsFixed(0),
      Value(value), Bounds(std::pair<double, double>(0, 0)),
      ErrType(ErrorType::NOTDEF), Error(std::pair<double, double>(0, 0)) {
  setError(error);
  setBounds(min, max);
}

DoubleParameter::DoubleParameter(const DoubleParameter &in)
    : Parameter(in.Name, ParType::DOUBLE) {
  *this = in;
}

void DoubleParameter::updateParameter(std::shared_ptr<DoubleParameter> newPar) {

  // Copy bounds
  if (newPar->hasBounds()) {
    HasBounds = 1;
    setBounds(newPar->bounds());
  } else
    HasBounds = 0;

  bool isFix = newPar->isFixed();
  fixParameter(0); // we ignore here if parameter is fixed

  // Copy value
  setValue(newPar->value());

  // Copy uncertainty
  if (newPar->errorType() == ErrorType::SYM)
    setError(newPar->error().first);
  else if (newPar->errorType() == ErrorType::ASYM)
    setError(newPar->error());
  else
    ErrType = ErrorType::NOTDEF;

  // Copy fix parameter
  fixParameter(isFix);
  return;
}

void DoubleParameter::setValue(const double inVal) {
  if (IsFixed)
    throw ParameterFixed("DoubleParameter:: () | Parameter " + name() +
                         " is fixed!");
  // Call notify only if value has changed! Otherwise tree is
  // recalculated also in case where current parameter is not changed
  if (Value == inVal)
    return;

  if (HasBounds && (inVal < bounds().first || inVal > bounds().second))
    throw ParameterOutOfBound("DoubleParameter::setValue() | Parameter " +
                              name() + " not within bounds: val=" +
                              std::to_string(inVal) + " [" +
                              std::to_string(bounds().first) + ";" +
                              std::to_string(bounds().second) + "]");

  Value = inVal;
  Notify();
}

std::pair<double, double> DoubleParameter::bounds() const { return Bounds; }

void DoubleParameter::setBounds(const double min, const double max) {
  if (!check_bounds(std::pair<double, double>(min, max)))
    throw BadParameter("DoubleParameter::setBounds() | Bounds no valid!");
  Bounds.first = min;
  Bounds.second = max;
}

void DoubleParameter::setBounds(const std::pair<double, double> r) {
  if (!check_bounds(r))
    throw BadParameter("DoubleParameter::setBounds() | Bounds no valid!");
  Bounds = r;
}

bool DoubleParameter::hasError() const {
  if (ErrType != ErrorType::NOTDEF)
    return true;
  return false;
}

std::pair<double, double> DoubleParameter::error() const {
  if (!hasError())
    throw std::runtime_error("DoubleParameter::error() | "
                             "Parameter " +
                             Name + " has no errors defined!");
  return Error;
}

void DoubleParameter::setError(double errLow, double errHigh) {
  SetErrorType(ErrorType::ASYM);
  Error.first = errLow;
  Error.second = errHigh;
}

void DoubleParameter::setError(std::pair<double, double> err) {
  SetErrorType(ErrorType::ASYM);
  Error = err;
}

void DoubleParameter::setError(double err) {
  SetErrorType(ErrorType::SYM);
  Error.first = err;
  Error.second = err;
}

bool DoubleParameter::operator==(const DoubleParameter otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->name() != otherPar.name())
    return false;
  if (this->value() != otherPar.value())
    return false;

  // We assume that if name and value are the same both parameters match. In
  // case that other properties of the parameterts differ we throw an
  // exception since we assume that this is a user mistake.
  if (HasBounds != otherPar.HasBounds || this->bounds() != otherPar.bounds())
    throw std::runtime_error("DoubleParameter::operator==() | Parameters "
                             "match by name (" +
                             name() + ") and value (" +
                             std::to_string(value()) +
                             ") but differs in "
                             "parameter bounds. We assume that there is a "
                             "mistake. Check your input files!");

  if (IsFixed != otherPar.IsFixed)
    throw std::runtime_error(
        "DoubleParameter::operator==() | Parameters "
        "match by name (" +
        name() + ") and value (" + std::to_string(value()) +
        ") but one is fixed the other not. "
        "We assume that there is a mistake. Check your input files!");

  if (ErrType != otherPar.ErrType || Error != otherPar.Error)
    throw std::runtime_error("DoubleParameter::operator==() | Parameters "
                             "match by name (" +
                             name() + ") and value (" +
                             std::to_string(value()) +
                             ") but differs in "
                             "parameter error. We assume that there is a "
                             "mistake. Check your input files!");

  return true;
}

bool DoubleParameter::check_bounds(
    const std::pair<double, double> bounds) const {
  if ((bounds.second > bounds.first) && (bounds.second >= Value) &&
      (bounds.first <= Value))
    return true;
  return false;
}

std::string DoubleParameter::to_str() const {
  std::stringstream oss;
  oss << Name;
  oss << "\t Val = " << Value;
  if (ErrType == ErrorType::SYM) {
    oss << " (+-" << Error.first << ")";
  } else if (ErrType == ErrorType::ASYM) {
    oss << " (+" << Error.second << " -" << Error.first << ")";
  }

  if (HasBounds)
    oss << "\t  [" << Bounds.first << " ; " << Bounds.second << "]";
  oss << " fix? " << isFixed();
  oss << "\t " << className();
  return oss.str();
}

std::string DoubleParameter::val_to_str() const {
  std::stringstream ovs;
  ovs << Value;
  return ovs.str();
}

void DoubleParameter::load(const boost::property_tree::ptree pt) {

  // Require that name and value are provided
  this->setName(pt.get<std::string>("<xmlattr>.Name"));
  this->setValue(pt.get<double>("Value"));

  // Optional settings
  if (pt.get_optional<double>("Error")) {
    this->setError(pt.get<double>("Error"));
  }
  if (pt.get_optional<double>("ErrorLow")) {
    if (pt.get_optional<double>("ErrorHigh"))
      this->setError(std::pair<double, double>(pt.get<double>("ErrorLow"),
                                               pt.get<double>("ErrorHigh")));
    else
      throw std::runtime_error("DoubleParameterFactory() | Parameter asymmetic "
                               "error not properly set!");
  } else if (pt.get_optional<double>("ErrorHigh")) {
    throw std::runtime_error("DoubleParameterFactory() | Parameter asymmetic "
                             "error not properly set!");
  } else { // Do not set a asymmetric errors
  }

  auto fix = pt.get_optional<bool>("Fix");
  if (fix)
    this->fixParameter(fix.get());
  else
    this->fixParameter(true);

  auto min = pt.get_optional<double>("Min");
  auto max = pt.get_optional<double>("Max");
  if (min && max) {
    this->setBounds(min.get(), max.get());
  } else if (min || max) { // Bounds not completely specified
    throw std::runtime_error(
        "DoubleParameterFactory() | Parameter bounds not properly set!");
  } else { // No bounds specified
  }
}

boost::property_tree::ptree DoubleParameter::save() const {
  boost::property_tree::ptree pt;

  // Require that name and value are provided
  pt.put("<xmlattr>.Name", name());
  pt.put("Value", value());
  pt.put("Fix", isFixed());
  if (hasBounds()) {
    pt.put("Min", bounds().first);
    pt.put("Max", bounds().second);
  }
  if (ErrType == ErrorType::SYM) {
    pt.put("Error", error().first);
  } else if (ErrType == ErrorType::ASYM) {
    pt.put("ErrorLow", error().first);
    pt.put("ErrorHigh", error().second);
  }
  return pt;
}

//================================= INTEGER ==================================
IntegerParameter::IntegerParameter(std::string inName)
    : Parameter(inName, ParType::INTEGER), Value(0), IsFixed(true),
      HasBounds(false), Bounds(std::pair<double, double>(0, 0)),
      ErrType(ErrorType::NOTDEF), Error(std::pair<double, double>(0, 0)) {}

IntegerParameter::IntegerParameter(std::string inName, const int value)
    : Parameter(inName, ParType::INTEGER), Value(value), IsFixed(true),
      HasBounds(false), Bounds(std::pair<double, double>(0, 0)),
      ErrType(ErrorType::NOTDEF), Error(std::pair<double, double>(0, 0)) {}

IntegerParameter::IntegerParameter(std::string inName, const int value,
                                   const int error)
    : Parameter(inName, ParType::INTEGER), Value(value), IsFixed(true),
      HasBounds(false), Bounds(std::pair<double, double>(0, 0)),
      ErrType(ErrorType::SYM), Error(std::pair<double, double>(error, error)) {}

IntegerParameter::IntegerParameter(std::string inName, const int value,
                                   const int min, const int max)
    : Parameter(inName, ParType::INTEGER), Value(value), IsFixed(true),
      HasBounds(false), Bounds(std::pair<double, double>(0, 0)),
      ErrType(ErrorType::SYM), Error(std::pair<double, double>(0, 0)) {

  if (check_bounds(min, max)) {
    Bounds = std::pair<double, double>(min, max);
    HasBounds = true;
  }
}

IntegerParameter::IntegerParameter(std::string inName, const int value,
                                   const int min, const int max,
                                   const int error)
    : Parameter(inName, ParType::INTEGER), Value(value), IsFixed(true),
      Bounds(std::pair<double, double>(0, 0)), ErrType(ErrorType::SYM),
      Error(std::pair<double, double>(error, error)) {

  setBounds(min, max);
}

IntegerParameter::IntegerParameter(const IntegerParameter &in)
    : Parameter(in.Name, ParType::INTEGER) {
  *this = in;
}

void IntegerParameter::setValue(const int inVal) {
  if (IsFixed) {
    throw ParameterFixed();
    return;
  }
  if (Value == inVal)
    return;
  Value = inVal;
  Notify();
}

bool IntegerParameter::hasError() const {
  if (ErrType != ErrorType::NOTDEF)
    return true;
  return false;
}

void IntegerParameter::setError(const int inErr) {
  ErrType = ErrorType::SYM;
  Error.first = inErr;
  Error.second = inErr;
}

void IntegerParameter::setError(const std::pair<int, int> inErr) {
  ErrType = ErrorType::ASYM;
  Error = inErr;
}

void IntegerParameter::setBounds(const int min, const int max) {
  if (check_bounds(min, max)) {
    Bounds = std::pair<double, double>(min, max);
    HasBounds = true;
  }
}

void IntegerParameter::setBounds(const std::pair<int, int> r) {
  if (check_bounds(r.first, r.second)) {
    Bounds = r;
    HasBounds = true;
  }
}

bool IntegerParameter::operator==(const IntegerParameter otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->name() != otherPar.name())
    return false;
  if (this->value() != otherPar.value())
    return false;

  // We assume that if name and value are the same both parameters match. In
  // case that other properties of the parameterts differ we throw an
  // exception since we assume that this is a user mistake.
  if (HasBounds != otherPar.HasBounds || Bounds != otherPar.Bounds)
    throw std::runtime_error("IntegerParameter::operator==() | Parameters "
                             "match by name (" +
                             name() + ") and value (" +
                             std::to_string(value()) +
                             ") but differs in "
                             "parameter bounds. We assume that there is a "
                             "mistake. Check your input files!");

  if (IsFixed != otherPar.IsFixed)
    throw std::runtime_error(
        "IntegerParameter::operator==() | Parameters "
        "match by name (" +
        name() + ") and value (" + std::to_string(value()) +
        ") but one is fixed the other not. "
        "We assume that there is a mistake. Check your input files!");

  if (Error != otherPar.Error || ErrType != otherPar.ErrType)
    throw std::runtime_error("IntegerParameter::operator==() | Parameters "
                             "match by name (" +
                             name() + ") and value (" +
                             std::to_string(value()) +
                             ") but differs in "
                             "parameter error. We assume that there is a "
                             "mistake. Check your input files!");

  return true;
}

bool IntegerParameter::check_bounds(const int min, const int max) const {
  if ((max > min) && (max >= Value) && (min <= Value))
    return true;
  return false;
}

std::string IntegerParameter::to_str() const {
  std::stringstream oss;
  oss << Name;
  oss << "\t Val = " << Value;
  oss << "\t Type = " << className();
  return oss.str();
}

std::string IntegerParameter::val_to_str() const {
  std::stringstream ovs;
  ovs << Value;
  return ovs.str();
}
//================================== BOOL ====================================

BoolParameter::BoolParameter(std::string inName)
    : Parameter(inName, ParType::BOOL), Value(0) {}

BoolParameter::BoolParameter(std::string inName, const bool value)
    : Parameter(inName, ParType::BOOL), Value(value) {}

BoolParameter::BoolParameter(const BoolParameter &in)
    : Parameter(in.Name, ParType::BOOL) {
  *this = in;
}

void BoolParameter::setValue(const bool inVal) {
  if (Value == inVal)
    return;
  Value = inVal;
  Notify();
}

bool BoolParameter::operator==(const BoolParameter otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->name() != otherPar.name())
    return false;
  if (this->value() != otherPar.value())
    return false;

  return true;
}

std::string BoolParameter::to_str() const {
  std::stringstream oss;
  oss << Name;
  oss << "\t Val = " << Value;
  oss << "\t Type = " << className();
  return oss.str();
}

std::string BoolParameter::val_to_str() const {
  std::stringstream ovs;
  ovs << Value;
  return ovs.str();
}
