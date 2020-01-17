// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// ComPWA exceptions.
///

#ifndef PWAEXCEPTIONS_H
#define PWAEXCEPTIONS_H

#include <exception>
#include <string>

namespace ComPWA {

///
/// \class Exception
/// ComPWA Exceptions base class.
/// This class defines the ComPWA exception base-class and provides a set of
/// standard exceptions.
///
class Exception : public std::exception {
public:
  Exception(const Exception &e) throw() : std::exception(e), what_(e.what_) {}

  Exception &operator=(const Exception &rhs) throw() {
    what_ = rhs.what_;
    return *this;
  }

  virtual ~Exception() throw() {}

  virtual const char *what() const throw() { return what_.c_str(); }

protected:
  Exception(const char *w = "") throw() : what_(w) {}

  Exception(const std::string &w) throw() : what_(w) {}

  std::string what_;
};

//------------------------------------------------------------------------------
//! @class   BadConfig
//!
//! @brief   Config is not complete
//------------------------------------------------------------------------------
class BadConfig : public Exception {
public:
  BadConfig(const std::string &error) : Exception(error) {}
  BadConfig(const char *error) : Exception(error) {}
  virtual ~BadConfig() throw() {}
};

//------------------------------------------------------------------------------
//! @class   BadParameter
//!
//! @brief   Parameter not existing
//------------------------------------------------------------------------------
class BadParameter : public Exception {
public:
  BadParameter(const std::string &error) : Exception(error) {}
  BadParameter(const char *error) : Exception(error) {}
  virtual ~BadParameter() throw() {}
};

//------------------------------------------------------------------------------
//! @class   BadIndex
//!
//! @brief   Index out of range
//------------------------------------------------------------------------------
class BadIndex : public Exception {
public:
  BadIndex(const std::string &error) : Exception(error) {}
  BadIndex(const char *error) : Exception(error) {}
  virtual ~BadIndex() throw() {}
};

//------------------------------------------------------------------------------
//! @class   CorruptFile
//!
//! @brief   Input data file is corrupt or incomplete
//------------------------------------------------------------------------------
class CorruptFile : public Exception {
public:
  CorruptFile(const std::string &error) : Exception(error) {}
  CorruptFile(const char *error) : Exception(error) {}
  virtual ~CorruptFile() throw() {}
};

//------------------------------------------------------------------------------
//! @class   ParameterFixed
//!
//! @brief   Parameter cannot be changed
//------------------------------------------------------------------------------
class ParameterFixed : public Exception {
public:
  ParameterFixed(const std::string &error = "Parameter fixed")
      : Exception(error) {}
  ParameterFixed(const char *error) : Exception(error) {}
  virtual ~ParameterFixed() throw() {}
};

//------------------------------------------------------------------------------
//! @class   WrongParType
//!
//! @brief   Parameter of wrong type
//------------------------------------------------------------------------------
class WrongParType : public Exception {
public:
  WrongParType(const std::string &error = "Parameter type wrong!")
      : Exception(error) {}
  WrongParType(const char *error) : Exception(error) {}
  virtual ~WrongParType() throw() {}
};

//------------------------------------------------------------------------------
//! @class   BeyondPhsp
//!
//! @brief   Data beyond phasespace requested
//------------------------------------------------------------------------------
class BeyondPhsp : public Exception {
public:
  BeyondPhsp(const std::string &error = "Data beyond phsp!")
      : Exception(error) {}
  BeyondPhsp(const char *error) : Exception(error) {}
  virtual ~BeyondPhsp() throw() {}
};

//------------------------------------------------------------------------------
//! @class   WrongVariableID
//!
//! @brief   Variable not found
//------------------------------------------------------------------------------
class WrongVariableID : public Exception {
public:
  WrongVariableID(const std::string &error = "Variable does not exist!")
      : Exception(error) {}
  WrongVariableID(const char *error) : Exception(error) {}
  virtual ~WrongVariableID() throw() {}
};
//------------------------------------------------------------------------------
//! @class   ParameterOutOfBound
//!
//! @brief   Parameter out of bound
//------------------------------------------------------------------------------
class ParameterOutOfBound : public Exception {
public:
  ParameterOutOfBound(
      const std::string &error = "Variable not within its limits!")
      : Exception(error) {}
  ParameterOutOfBound(const char *error) : Exception(error) {}
  virtual ~ParameterOutOfBound() throw() {}
};
//------------------------------------------------------------------------------
//! @class   TreeBuildError
//!
//! @brief   Error in tree at construction
//------------------------------------------------------------------------------
class TreeBuildError : public Exception {
public:
  TreeBuildError(const std::string &error = "Tree can not be build!")
      : Exception(error) {}
  TreeBuildError(const char *error) : Exception(error) {}
  virtual ~TreeBuildError() throw() {}
};

} /* namespace ComPWA */
#endif
//******************************************************************************
//! EOF
//******************************************************************************
