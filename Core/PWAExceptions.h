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
//! @brief   Cannot connect to database or config is not complete
//------------------------------------------------------------------------------
class BadConfig : public Exception {
public:
  BadConfig ( const std::string& error ) :
    Exception(error)
  {}
  BadConfig ( const char *error ) :
    Exception(error)
  {}
  ~BadConfig () throw() {}
};


#endif

//******************************************************************************
//! EOF
//******************************************************************************
