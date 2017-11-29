// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Parameter.hpp"

using namespace ComPWA;

//=============================== MULTICOMPEX ================================
//! Getter for value of parameter
std::complex<double> MultiComplex::GetValue(unsigned int i) const {
  if (i >= val_.size())
    return 0;
  return val_[i];
}

//! Setter for value of parameter
void MultiComplex::SetValue(const std::complex<double> inVal, unsigned int i) {
  if (i >= val_.size())
    return;
  if (val_[i] == inVal)
    return;
  val_[i] = inVal;
  Notify();
}

bool MultiComplex::operator==(const MultiComplex otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->GetName() != otherPar.GetName())
    return false;
  if (this->GetValue() != otherPar.GetValue())
    return false;

  return true;
}

std::string MultiComplex::make_str() const {
  std::stringstream oss;
  oss << name_;
  size_t max = val_.size();
  if (max > 5)
    max = 5; // display only 5 variables
  oss << "\t Val = ";
  for (unsigned int i = 0; i < max - 1; i++)
    oss << val_[i] << ", ";
  oss << val_[max - 1];
  if (max < val_.size())
    oss << " ... ";
  oss << "\t Type = " << TypeName();
  return oss.str();
}

std::string MultiComplex::make_val_str() const {
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

double MultiDouble::GetValue(unsigned int i) const {
  if (i >= val_.size())
    return 0;
  return val_[i];
}

void MultiDouble::SetValue(const double inVal, unsigned int i) {
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
  if (this->GetName() != otherPar.GetName())
    return false;
  if (this->GetValue() != otherPar.GetValue())
    return false;

  return true;
}

std::string MultiDouble::make_str() const {
  std::stringstream oss;
  oss << name_;
  size_t max = val_.size();
  if (max > 5)
    max = 5; // display only 10 variables
  oss << "\t Val = ";
  for (unsigned int i = 0; i < max - 1; i++)
    oss << val_[i] << ", ";
  oss << val_[max - 1];
  if (max < val_.size())
    oss << " ... ";
  oss << "\t Type = " << TypeName();
  return oss.str();
}

std::string MultiDouble::make_val_str() const {
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
unsigned int MultiUnsignedInteger::GetValue(unsigned int i) const {
  if (i >= val_.size())
    return 0;
  return val_[i];
}

//! Setter for value of parameter
void MultiUnsignedInteger::SetValue(const unsigned int inVal, unsigned int i) {
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
  if (this->GetName() != otherPar.GetName())
    return false;
  if (this->GetValue() != otherPar.GetValue())
    return false;

  return true;
}
std::string MultiUnsignedInteger::make_str() const {
  std::stringstream oss;
  oss << name_;
  size_t max = val_.size();
  if (max > 5)
    max = 5; // display only 10 variables
  oss << "\t Val = ";
  for (unsigned int i = 0; i < max - 1; i++)
    oss << val_[i] << ", ";
  oss << val_[max - 1];
  if (max < val_.size())
    oss << " ... ";
  oss << "\t Type = " << TypeName();
  return oss.str();
}

std::string MultiUnsignedInteger::make_val_str() const {
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
    : Parameter(inName, ParType::COMPLEX), val_(0., 0.), err_(0., 0.),
      min_(0., 0.), max_(0., 0.) {
  bounds_ = usebounds_ = hasError_ = fixed_ = false;
}

ComplexParameter::ComplexParameter(std::string inName,
                                   const std::complex<double> value)
    : Parameter(inName, ParType::COMPLEX), val_(value), err_(0, 0),
      min_(0, 0), max_(0, 0) {
  bounds_ = usebounds_ = hasError_ = fixed_ = false;
}

ComplexParameter::ComplexParameter(std::string inName,
                                   const std::complex<double> value,
                                   const std::complex<double> error)
    : Parameter(inName, ParType::COMPLEX), val_(value), err_(error),
      min_(0, 0), max_(0, 0) {
  bounds_ = usebounds_ = fixed_ = false;
  hasError_ = true;
}

ComplexParameter::ComplexParameter(std::string inName,
                                   const std::complex<double> value,
                                   const std::complex<double> min,
                                   const std::complex<double> max)
    : Parameter(inName, ParType::COMPLEX), val_(value), err_(0, 0),
      min_(0, 0), max_(0, 0) {
  bounds_ = usebounds_ = hasError_ = fixed_ = false;
  if (check_bounds(min, max)) {
    min_ = min;
    max_ = max;
    bounds_ = true;
  }
}
ComplexParameter::ComplexParameter(std::string inName,
                                   const std::complex<double> value,
                                   const std::complex<double> min,
                                   const std::complex<double> max,
                                   const std::complex<double> error)
    : Parameter(inName, ParType::COMPLEX), val_(value), err_(error),
      min_(0, 0), max_(0, 0) {
  bounds_ = usebounds_ = fixed_ = false;
  hasError_ = true;
  if (check_bounds(min, max)) {
    min_ = min;
    max_ = max;
    bounds_ = true;
  }
}

ComplexParameter::ComplexParameter(const ComplexParameter &in)
    : Parameter(in.name_, ParType::COMPLEX) {
  *this = in;
}

void ComplexParameter::SetValue(const std::complex<double> inVal) {
  if (fixed_) {
    throw ParameterFixed();
    return;
  }
  if (val_ == inVal)
    return;
  val_ = inVal;
  Notify();
}

void ComplexParameter::SetError(const std::complex<double> inErr) {
  err_ = inErr;
  hasError_ = true;
}

bool ComplexParameter::SetMinMax(const std::complex<double> inMin,
                                 const std::complex<double> inMax) {
  bool valid = check_bounds(inMin, inMax);
  if (valid) {
    min_ = inMin;
    max_ = inMax;
    bounds_ = true;
  }
  return valid;
}

bool ComplexParameter::SetMinValue(const std::complex<double> min) {
  bool valid = check_bounds(min, max_);
  if (valid) {
    min_ = min;
    bounds_ = true;
  }
  return valid;
}

bool ComplexParameter::SetMaxValue(const std::complex<double> max) {
  bool valid = check_bounds(min_, max);
  if (valid) {
    max_ = max;
    bounds_ = true;
  }
  return valid;
}
bool ComplexParameter::operator==(const ComplexParameter &otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->GetName() != otherPar.GetName())
    return false;
  if (this->GetValue() != otherPar.GetValue())
    return false;

  /* We assume that if name and value are the same both parameters match. In
   * case that other properties of the parameterts differ we throw an
   * exception since we assume that this is a user mistake.
   */
  if (bounds_ != otherPar.bounds_ || usebounds_ != otherPar.usebounds_ ||
      this->GetMinValue() != otherPar.GetMinValue() ||
      this->GetMaxValue() != otherPar.GetMaxValue())
    throw std::runtime_error("ComplexParameter::operator==() | Parameters "
                             "match by name (" +
                             GetName() + ") and value (" +
                             std::to_string(GetValue().real()) + "+i*" +
                             std::to_string(GetValue().imag()) +
                             ") but differs in "
                             "parameter bounds. We assume that there is a "
                             "mistake. Check your input files!");

  if (fixed_ != otherPar.fixed_)
    throw std::runtime_error(
        "ComplexParameter::operator==() | Parameters "
        "match by name (" +
        GetName() + ") and value (" + std::to_string(GetValue().real()) +
        "+i*" + std::to_string(GetValue().imag()) +
        ") but one is fixed the other not. "
        "We assume that there is a mistake. Check your input files!");

  if (hasError_ != otherPar.hasError_ || err_ != otherPar.err_)
    throw std::runtime_error("ComplexParameter::operator==() | Parameters "
                             "match by name (" +
                             GetName() + ") and value (" +
                             std::to_string(GetValue().real()) + "+i*" +
                             std::to_string(GetValue().imag()) +
                             ") but differs in "
                             "parameter error. We assume that there is a "
                             "mistake. Check your input files!");

  return true;
}

bool ComplexParameter::check_bounds(const std::complex<double> min,
                                    const std::complex<double> max) {
  if ((max.real() > min.real()) && (max.real() >= val_.real()) &&
      (min.real() <= val_.real()) && (max.imag() > min.imag()) &&
      (max.imag() >= val_.imag()) && (min.imag() <= val_.imag()))
    return true;
  return false;
}

std::string ComplexParameter::make_str() const {
  std::stringstream oss;
  oss << name_;
  oss << "\t Val = " << val_;
  if (bounds_)
    oss << "\t  Min-Max = " << min_ << " to " << max_;
  if (hasError_)
    oss << "\t  Err = " << err_;
  oss << "\t Type = " << TypeName();
  return oss.str();
}

std::string ComplexParameter::make_val_str() const {
  std::stringstream ovs;
  ovs << val_;
  return ovs.str();
}
//================================= DOUBLE ===================================

DoubleParameter::DoubleParameter(std::string inName)
    : Parameter(inName, ParType::DOUBLE), bounds_(false), fixed_(0), val_(0),
      min_(0), max_(0), errorType(ErrorType::NOTDEF), errorLow(0.),
      errorHigh(0.) {}

DoubleParameter::DoubleParameter(std::string inName, const double value)
    : Parameter(inName, ParType::DOUBLE), bounds_(false), fixed_(0),
      val_(value), min_(0), max_(0), errorType(ErrorType::NOTDEF), errorLow(0.),
      errorHigh(0.) {}

DoubleParameter::DoubleParameter(std::string inName, const double value,
                                 const double error)
    : Parameter(inName, ParType::DOUBLE), bounds_(false), fixed_(0),
      val_(value), min_(0), max_(0), errorType(ErrorType::NOTDEF), errorLow(0.),
      errorHigh(0.) {
  SetError(error);
}

DoubleParameter::DoubleParameter(std::string inName, const double value,
                                 const double min, const double max)
    : Parameter(inName, ParType::DOUBLE), bounds_(false), fixed_(0),
      val_(value), min_(0), max_(0), errorType(ErrorType::NOTDEF), errorLow(0.),
      errorHigh(0.) {
  SetMinMax(min, max);
}

DoubleParameter::DoubleParameter(std::string inName, const double value,
                                 const double min, const double max,
                                 const double error)
    : Parameter(inName, ParType::DOUBLE), bounds_(false), fixed_(0),
      val_(value), min_(0), max_(0), errorType(ErrorType::NOTDEF), errorLow(0.),
      errorHigh(0.) {
  SetError(error);
  SetMinMax(min, max);
}

DoubleParameter::DoubleParameter(const DoubleParameter &in)
    : Parameter(in.name_, ParType::DOUBLE) {
  *this = in;
}

void DoubleParameter::UpdateParameter(std::shared_ptr<DoubleParameter> newPar) {

  // Copy bounds
  if (newPar->HasBounds()){
    bounds_ = 1;
    SetMinMax(newPar->GetMinValue(), newPar->GetMaxValue());
  } else
    bounds_ = 0;

  bool isFix = newPar->IsFixed();
  FixParameter(0); // we ignore here if parameter is fixed
  
  // Copy value
  SetValue(newPar->GetValue());

  // Copy uncertainty
  if (newPar->GetErrorType() == ErrorType::SYM)
    SetError(newPar->GetError());
  else if (newPar->GetErrorType() == ErrorType::ASYM)
    SetError(newPar->GetErrorLow(), newPar->GetErrorHigh());
  else
    SetErrorType(ErrorType::NOTDEF);

  // Copy fix parameter
  FixParameter(isFix);
  return;
}

void DoubleParameter::SetValue(const double inVal) {
  if (fixed_)
    throw ParameterFixed("DoubleParameter::SetValue() | Parameter " +
                         GetName() + " is fixed!");
  // Call notify only if value has changed! Otherwise tree is
  // recalculated also in case where current parameter is not changed
  if (val_ == inVal)
    return;

  if (bounds_ && (inVal < GetMinValue() || inVal > GetMaxValue()))
    throw ParameterOutOfBound(
        "DoubleParameter::SetValue() | "
        "Parameter " +
        GetName() + " not within bounds: val=" + std::to_string(inVal) + " [" +
        std::to_string((long double)GetMinValue()) + ";" +
        std::to_string((long double)GetMaxValue()) + "]!");

  val_ = inVal;
  Notify();
}

void DoubleParameter::SetRange(const double min, const double max) {
  SetMinMax(min, max);
}

void DoubleParameter::SetMinMax(const double min, const double max) {
  try {
    SetMinValue(min);
  } catch (ParameterOutOfBound &ex) {
  }
  try {
    SetMaxValue(max);
  } catch (ParameterOutOfBound &ex) {
  }
  if (!check_bounds(min_, max_))
    throw ParameterOutOfBound("DoubleParameter::SetMinMaxValue() | "
                             "Bounds not valid for parameter " +
                             GetName() + ": " + std::to_string(GetValue()) +
                             " [" + std::to_string((long double)min_) + ";" +
                             std::to_string((long double)max_) + "]!");
  bounds_ = true;
}

void DoubleParameter::SetMinValue(const double min) {
  min_ = min;
  if (!check_bounds(min_, max_))
    throw ParameterOutOfBound("DoubleParameter::SetMinValue() | "
                              "Boundary not valid: [" +
                              std::to_string(GetMinValue()) + ", " +
                              std::to_string(GetMaxValue()) + "]!");
  bounds_ = true;
}
void DoubleParameter::SetMaxValue(const double max) {
  max_ = max;
  if (!check_bounds(min_, max_))
    throw ParameterOutOfBound("DoubleParameter::SetMaxValue() | "
                              "Boundary not valid: [" +
                              std::to_string(GetMinValue()) + ", " +
                              std::to_string(GetMaxValue()) + "]!");
  bounds_ = true;
}

bool DoubleParameter::HasError() const {
  if (GetErrorType() == ErrorType::NOTDEF)
    return 0;
  else
    return 1;
}

double DoubleParameter::GetError() const {
  if (!HasError())
    throw std::runtime_error("DoubleParameter::GetError() | "
                             "Parameter " +
                             name_ + " has no errors defined!");
  return (GetErrorHigh() + GetErrorLow()) / 2;
}

double DoubleParameter::GetErrorHigh() const {
  if (!HasError())
    throw std::runtime_error("DoubleParameter::GetError() | "
                             "Parameter " +
                             name_ + " has no errors defined!");
  //		if(GetErrorType()==ErrorType::SYM){
  //			LOG(info) << "DoubleParameter::GetErrorHigh() |
  // Parameter
  //"<<name_
  //					<<" has no asymmetric errors! Returning
  // symmetric
  // error";
  //			return GetError();
  //		}
  return errorHigh;
}

double DoubleParameter::GetErrorLow() const {
  GetName();
  if (!HasError())
    throw std::runtime_error("DoubleParameter::GetError() | "
                             "Parameter " +
                             name_ + " has no errors defined!");
  //		if(GetErrorType()==ErrorType::SYM){
  //			LOG(info) << "DoubleParameter::GetErrorHigh() |
  // Parameter
  //"<<name_
  //					<<" has no assymetric errors! Returning
  // symmetric
  // error";
  //			return GetError();
  //		}
  if (!HasError())
    throw std::runtime_error("DoubleParameter::GetError() | "
                             "Parameter " +
                             name_ + " has no errors defined!");
  return errorLow;
}

void DoubleParameter::SetError(double errLow, double errHigh) {
  // if(fixed_)
  //	throw ParameterFixed("DoubleParameter::SetError(double, double) |
  // Parameter "+GetName()+" is fixed!");
  SetErrorType(ErrorType::ASYM);
  SetErrorHigh(errHigh);
  SetErrorLow(errLow);
}

void DoubleParameter::SetError(double err) {
  // if(fixed_)
  //	throw ParameterFixed("DoubleParameter::SetError(double) | Parameter
  //"+GetName()+" is fixed!");
  SetErrorType(ErrorType::SYM);
  SetErrorHigh(err);
  SetErrorLow(err);
}

bool DoubleParameter::operator==(const DoubleParameter otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->GetName() != otherPar.GetName())
    return false;
  if (this->GetValue() != otherPar.GetValue())
    return false;

  /* We assume that if name and value are the same both parameters match. In
   * case that other properties of the parameterts differ we throw an
   * exception since we assume that this is a user mistake.
   */
  if (bounds_ != otherPar.bounds_ ||
      this->GetMinValue() != otherPar.GetMinValue() ||
      this->GetMaxValue() != otherPar.GetMaxValue())
    throw std::runtime_error("DoubleParameter::operator==() | Parameters "
                             "match by name (" +
                             GetName() + ") and value (" +
                             std::to_string(GetValue()) +
                             ") but differs in "
                             "parameter bounds. We assume that there is a "
                             "mistake. Check your input files!");

  if (fixed_ != otherPar.fixed_)
    throw std::runtime_error(
        "DoubleParameter::operator==() | Parameters "
        "match by name (" +
        GetName() + ") and value (" + std::to_string(GetValue()) +
        ") but one is fixed the other not. "
        "We assume that there is a mistake. Check your input files!");

  if (errorType != otherPar.errorType || errorLow != otherPar.errorLow ||
      errorHigh != otherPar.errorHigh)
    throw std::runtime_error("DoubleParameter::operator==() | Parameters "
                             "match by name (" +
                             GetName() + ") and value (" +
                             std::to_string(GetValue()) +
                             ") but differs in "
                             "parameter error. We assume that there is a "
                             "mistake. Check your input files!");

  return true;
}

bool DoubleParameter::check_bounds(const double min, const double max) {
  if ((max > min) && (max >= val_) && (min <= val_))
    return true;
  return false;
}

std::string DoubleParameter::make_str() const {
  std::stringstream oss;
  oss << name_;
  oss << "\t Val = " << val_;
  if (errorLow != errorHigh)
    oss << " (+" << errorHigh << " -" << errorLow << ")";
  else if (errorLow != 0)
    oss << " (+-" << errorLow << ")";
  if (bounds_)
    oss << "\t  [" << min_ << " ; " << max_ << "]";
  oss << " fix? " << IsFixed();
  oss << "\t " << TypeName();
  return oss.str();
}

std::string DoubleParameter::make_val_str() const {
  std::stringstream ovs;
  ovs << val_;
  return ovs.str();
}

//================================= INTEGER ==================================
IntegerParameter::IntegerParameter(std::string inName)
    : Parameter(inName, ParType::INTEGER), val_(0), min_(0), max_(0),
      err_(0) {
  bounds_ = usebounds_ = hasError_ = fixed_ = false;
}

IntegerParameter::IntegerParameter(std::string inName, const int value)
    : Parameter(inName, ParType::INTEGER), val_(value), min_(0), max_(0),
      err_(0) {
  bounds_ = usebounds_ = hasError_ = fixed_ = false;
}

IntegerParameter::IntegerParameter(std::string inName, const int value,
                                   const int error)
    : Parameter(inName, ParType::INTEGER), val_(value), min_(0), max_(0),
      err_(error) {
  bounds_ = usebounds_ = fixed_ = false;
  hasError_ = true;
}

IntegerParameter::IntegerParameter(std::string inName, const int value,
                                   const int min, const int max)
    : Parameter(inName, ParType::INTEGER), val_(value), min_(0), max_(0),
      err_(0) {
  bounds_ = usebounds_ = hasError_ = fixed_ = false;
  if (check_bounds(min, max)) {
    min_ = min;
    max_ = max;
    bounds_ = true;
  }
}

IntegerParameter::IntegerParameter(std::string inName, const int value,
                                   const int min, const int max,
                                   const int error)
    : Parameter(inName, ParType::INTEGER), val_(value), min_(0), max_(0),
      err_(error) {
  bounds_ = usebounds_ = fixed_ = false;
  hasError_ = true;
  if (check_bounds(min, max)) {
    min_ = min;
    max_ = max;
    bounds_ = true;
  }
}

IntegerParameter::IntegerParameter(const IntegerParameter &in)
    : Parameter(in.name_, ParType::INTEGER) {
  *this = in;
}

void IntegerParameter::SetValue(const int inVal) {
  if (fixed_) {
    throw ParameterFixed();
    return;
  }
  if (val_ == inVal)
    return;
  val_ = inVal;
  Notify();
}

void IntegerParameter::SetError(const int inErr) {
  err_ = inErr;
  hasError_ = true;
}

bool IntegerParameter::SetMinMax(const int inMin, const int inMax) {
  bool valid = check_bounds(inMin, inMax);
  if (valid) {
    min_ = inMin;
    max_ = inMax;
    bounds_ = true;
  }
  return valid;
}

bool IntegerParameter::SetMinValue(const int min) {
  bool valid = check_bounds(min, max_);
  if (valid) {
    min_ = min;
    bounds_ = true;
  }
  return valid;
}

bool IntegerParameter::SetMaxValue(const int max) {
  bool valid = check_bounds(min_, max);
  if (valid) {
    max_ = max;
    bounds_ = true;
  }
  return valid;
}

bool IntegerParameter::operator==(const IntegerParameter otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->GetName() != otherPar.GetName())
    return false;
  if (this->GetValue() != otherPar.GetValue())
    return false;

  /* We assume that if name and value are the same both parameters match. In
   * case that other properties of the parameterts differ we throw an
   * exception since we assume that this is a user mistake.
   */
  if (bounds_ != otherPar.bounds_ || usebounds_ != otherPar.usebounds_ ||
      this->GetMinValue() != otherPar.GetMinValue() ||
      this->GetMaxValue() != otherPar.GetMaxValue())
    throw std::runtime_error("IntegerParameter::operator==() | Parameters "
                             "match by name (" +
                             GetName() + ") and value (" +
                             std::to_string(GetValue()) +
                             ") but differs in "
                             "parameter bounds. We assume that there is a "
                             "mistake. Check your input files!");

  if (fixed_ != otherPar.fixed_)
    throw std::runtime_error(
        "IntegerParameter::operator==() | Parameters "
        "match by name (" +
        GetName() + ") and value (" + std::to_string(GetValue()) +
        ") but one is fixed the other not. "
        "We assume that there is a mistake. Check your input files!");

  if (err_ != otherPar.err_)
    throw std::runtime_error("IntegerParameter::operator==() | Parameters "
                             "match by name (" +
                             GetName() + ") and value (" +
                             std::to_string(GetValue()) +
                             ") but differs in "
                             "parameter error. We assume that there is a "
                             "mistake. Check your input files!");

  return true;
}

bool IntegerParameter::check_bounds(const int min, const int max) {
  if ((max > min) && (max >= val_) && (min <= val_))
    return true;
  return false;
}

std::string IntegerParameter::make_str() const {
  std::stringstream oss;
  oss << name_;
  oss << "\t Val = " << val_;
  if (bounds_)
    oss << "\t  Min-Max = " << min_ << " to " << max_;
  if (hasError_)
    oss << "\t  Err = " << err_;
  oss << "\t Type = " << TypeName();
  return oss.str();
}

std::string IntegerParameter::make_val_str() const {
  std::stringstream ovs;
  ovs << val_;
  return ovs.str();
}
//================================== BOOL ====================================

BoolParameter::BoolParameter(std::string inName)
    : Parameter(inName, ParType::BOOL), val_(0), err_(0) {
  hasError_ = fixed_ = usebounds_ = false;
}

BoolParameter::BoolParameter(std::string inName, const bool value)
    : Parameter(inName, ParType::BOOL), val_(value), err_(0) {
  hasError_ = fixed_ = usebounds_ = false;
}

BoolParameter::BoolParameter(std::string inName, const bool value,
                             const bool error)
    : Parameter(inName, ParType::BOOL), val_(value), err_(error) {
  usebounds_ = false;
  fixed_ = false;
  hasError_ = true;
}

BoolParameter::BoolParameter(const BoolParameter &in)
    : Parameter(in.name_, ParType::BOOL) {
  *this = in;
}

void BoolParameter::SetValue(const bool inVal) {
  if (fixed_) {
    throw ParameterFixed();
    return;
  }
  if (val_ == inVal)
    return;
  val_ = inVal;
  Notify();
}

void BoolParameter::SetError(const bool inErr) {
  err_ = inErr;
  hasError_ = true;
}

bool BoolParameter::operator==(const BoolParameter otherPar) const {
  if (this->type() != otherPar.type())
    return false;
  if (this->GetName() != otherPar.GetName())
    return false;
  if (this->GetValue() != otherPar.GetValue())
    return false;

  if (fixed_ != otherPar.fixed_)
    throw std::runtime_error(
        "BoolParameter::operator==() | Parameters "
        "match by name (" +
        GetName() + ") and value (" + std::to_string(GetValue()) +
        ") but one is fixed the other not. "
        "We assume that there is a mistake. Check your input files!");

  if (err_ != otherPar.err_)
    throw std::runtime_error("BoolParameter::operator==() | Parameters "
                             "match by name (" +
                             GetName() + ") and value (" +
                             std::to_string(GetValue()) +
                             ") but differs in "
                             "parameter error. We assume that there is a "
                             "mistake. Check your input files!");

  return true;
}

std::string BoolParameter::make_str() const {
  std::stringstream oss;
  oss << name_;
  oss << "\t Val = " << val_;
  if (hasError_)
    oss << "\t  Err = " << err_;
  oss << "\t Type = " << TypeName();
  return oss.str();
}

std::string BoolParameter::make_val_str() const {
  std::stringstream ovs;
  ovs << val_;
  return ovs.str();
}
