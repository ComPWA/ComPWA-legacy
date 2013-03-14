//! ComPWA Exceptions
/*! \class PWAException
 * @file PWAExceptions.hpp
 * This class defines the ComPWA exception base-class and provides a set of
 * standard exceptions.
*/

#ifndef PWAEXCEPTIONS_H
#define PWAEXCEPTIONS_H

#include <exception>
#include <string>


class PWAException : public std::exception {
public:
  PWAException(const PWAException& e) throw() :
    std::exception(e),
    what_(e.what_)
  {}

  PWAException& operator=(const PWAException& rhs) throw() {
    what_ = rhs.what_;
    return *this;
  }

  virtual ~PWAException() throw() { }

  virtual const char* what() const throw() {
    return what_.c_str();
  }

protected:
  PWAException(const char* w = "") throw() :
    what_(w)
  {}

  PWAException(const std::string& w) throw() :
    what_(w)
  {}

  std::string what_;
};

//------------------------------------------------------------------------------
//! @class   BadConfig
//!
//! @brief   Config is not complete
//------------------------------------------------------------------------------
class BadConfig : public PWAException {
public:
  BadConfig ( const std::string& error ) :
    PWAException(error)
  {}
  BadConfig ( const char *error ) :
    PWAException(error)
  {}
  virtual ~BadConfig () throw() {}
};

//------------------------------------------------------------------------------
//! @class   BadParameter
//!
//! @brief   Parameter not existing
//------------------------------------------------------------------------------
class BadParameter : public PWAException {
public:
  BadParameter ( const std::string& error ) :
    PWAException(error)
  {}
  BadParameter ( const char *error ) :
    PWAException(error)
  {}
  virtual ~BadParameter () throw() {}
};

//------------------------------------------------------------------------------
//! @class   ParameterFixed
//!
//! @brief   Parameter cannot be changed
//------------------------------------------------------------------------------
class ParameterFixed : public PWAException {
public:
  ParameterFixed ( const std::string& error = "Parameter fixed" ) :
    PWAException(error)
  {}
  ParameterFixed ( const char *error ) :
    PWAException(error)
  {}
  virtual ~ParameterFixed () throw() {}
};


#endif

//******************************************************************************
//! EOF
//******************************************************************************
