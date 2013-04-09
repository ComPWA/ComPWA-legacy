//! ComPWA Exceptions
/*! \class Exception
 * @file PWAExceptions.hpp
 * This class defines the ComPWA exception base-class and provides a set of
 * standard exceptions.
*/

#ifndef PWAEXCEPTIONS_H
#define PWAEXCEPTIONS_H

#include <exception>
#include <string>


class Exception : public std::exception {
public:
  Exception(const Exception& e) throw() :
    std::exception(e),
    what_(e.what_)
  {}

  Exception& operator=(const Exception& rhs) throw() {
    what_ = rhs.what_;
    return *this;
  }

  virtual ~Exception() throw() { }

  virtual const char* what() const throw() {
    return what_.c_str();
  }

protected:
  Exception(const char* w = "") throw() :
    what_(w)
  {}

  Exception(const std::string& w) throw() :
    what_(w)
  {}

  std::string what_;
};

//------------------------------------------------------------------------------
//! @class   BadConfig
//!
//! @brief   Config is not complete
//------------------------------------------------------------------------------
class BadConfig : public Exception {
public:
  BadConfig ( const std::string& error ) :
    Exception(error)
  {}
  BadConfig ( const char *error ) :
    Exception(error)
  {}
  virtual ~BadConfig () throw() {}
};

//------------------------------------------------------------------------------
//! @class   BadParameter
//!
//! @brief   Parameter not existing
//------------------------------------------------------------------------------
class BadParameter : public Exception {
public:
  BadParameter ( const std::string& error ) :
    Exception(error)
  {}
  BadParameter ( const char *error ) :
    Exception(error)
  {}
  virtual ~BadParameter () throw() {}
};

//------------------------------------------------------------------------------
//! @class   ParameterFixed
//!
//! @brief   Parameter cannot be changed
//------------------------------------------------------------------------------
class ParameterFixed : public Exception {
public:
  ParameterFixed ( const std::string& error = "Parameter fixed" ) :
    Exception(error)
  {}
  ParameterFixed ( const char *error ) :
    Exception(error)
  {}
  virtual ~ParameterFixed () throw() {}
};


#endif

//******************************************************************************
//! EOF
//******************************************************************************
