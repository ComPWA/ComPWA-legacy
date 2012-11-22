//! Base class for internal parameter.
/*! \class PWAParameter
 * @file PWAParameter.hpp
 * This class defines the interface for the internal container of a fit
 * parameter. A parameter consists of a value with optional bounds and error.
*/

#ifndef _PWAPARAMETER_HPP_
#define _PWAPARAMETER_HPP_

#include <iostream>
#include <string>
#include <sstream>

class PWAParameter
{

public:

  //! Standard constructor without information
  /*!
   * Standard constructor without information provided.
  */
  PWAParameter() {
  }

  //! Empty Destructor
  /*!
   * There is nothing to destroy :(
  */
  virtual ~PWAParameter() { /* nothing */	}

  //! Check if parameter has bounds
  virtual const inline bool HasBounds() const {return bounds_;}
  //! Check if bounds should be used
  virtual const inline bool UseBounds() const {if(bounds_)return usebounds_; return false;}
  //! Check if parameter has an error
  virtual const inline bool HasError() const {return error_;}

  //! Getter for value of parameter cast to double
  virtual const double GetValue() const = 0;
  //! Getter for lower bound of parameter cast to double
  virtual const double GetMinValue() const = 0;
  //! Getter for upper bound of parameter cast to double
  virtual const double GetMaxValue() const = 0;
  //! Getter for error of parameter cast to double
  virtual const double GetError() const = 0;

  //! Setter for value of parameter cast to double
  virtual void SetValue(const double inVal) = 0;
  //! Setter for bounds of parameter cast to double
  virtual const bool SetMinMax(const double inMin, const double inMax) = 0;
  //! Setter for error of parameter cast to double
  virtual void SetError(const double inErr) = 0;

  //! Set if bounds should be used
  virtual const inline void UseBounds(const bool use) {usebounds_=use;}

  //! A public function returning a string with parameter information
  /*!
   * This function simply returns the member string out_, which contains
   * all parameter information. The string gets rebuild with every change
   * of the parameter.
   * \return string with parameter information
   * \sa operator<<, make_str()
  */
  std::string const& to_str() const{
    return out_;
  }

  //! A public function returning a string naming its type
  /*!
   * This function is used to get the type of the implementation of this
   * general parameter interface.
  */
  virtual const std::string type() = 0;

protected:
  std::string out_; /*!< Output string to print information */
  bool bounds_; /*!< Are valid bounds defined for this parameter? */
  bool error_; /*!< Is an error defined for this parameter? */
  bool usebounds_; /*!< Do you want to restrict your parameter? */

  //! A protected function which creates an output string for printing
  /*!
   * This function uses all available information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str()
  */
  virtual void make_str() = 0;

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type. The definition
   * of this class has to be outside the namespace of the class.
   * \sa make_str(), to_str()
  */
  friend std::ostream & operator<<(std::ostream &os, const PWAParameter& p);

};

inline std::ostream & operator<<(std::ostream &os, const PWAParameter& p){
  return os << p.to_str();
}

#endif
