//! Internal container representing a generic parameter.
/*! \class PWAGenericPar
 * @file PWAGenericPar.hpp
 * This class provides a internal container for information of a fit parameter.
 * A parameter consists of a value with optional bounds and error. It is a
 * template so it can represent integer, floating-point or boolean values.
 * It implements the PWAParameter interface and therefore represents parameters
 * for ComPWA with generic type's.
*/

#ifndef _PWAGENERICPAR_HPP_
#define _PWAGENERICPAR_HPP_

#include <iostream>
#include <string>
#include <sstream>
#include <typeinfo>

#include "Core/PWAParameter.hpp"

template <class T>
class PWAGenericPar : public PWAParameter
{

public:

  //! Standard constructor without information
  /*!
   * Standard constructor with no information provided. Creates parameter
   * with value 0 but without bounds or an error.
   * \sa make_str()
  */
  PWAGenericPar():val_(0),min_(0),max_(0),err_(0) {
    bounds_= usebounds_ = error_ = false;
    make_str();
  }

  //! Standard constructor with a value
  /*!
   * Standard constructor with just a value provided. Creates parameter
   * with given value but without bounds or an error.
   * \param value input value of the parameter
   * \sa make_str()
  */
  PWAGenericPar(const T value):val_(value),min_(0),max_(0),err_(0){
    bounds_= usebounds_ = error_ = false;
    make_str();
  }

  //! Standard constructor with value and error
  /*!
   * Standard constructor with value and error provided. Creates parameter
   * with given value and error but without bounds.
   * \param value input value of the parameter
   * \param error input error of the parameter
   * \sa make_str()
  */
  PWAGenericPar(const T value, const T error)
  :val_(value),min_(0),max_(0),err_(error){
    bounds_= usebounds_ = false;
    error_ = true;
    make_str();
  }

  //! Standard constructor with value and bounds
  /*!
   * Standard constructor with value and bounds provided. Creates parameter
   * with given value and bounds but without error. If a check for valid
   * bounds fails, just the value is used.
   * \param value input value of the parameter
   * \param min input lower bound
   * \param max input upper bound
   * \sa make_str(), check_bounds()
  */
  PWAGenericPar(const T value, const T min, const T max)
  :val_(value),min_(0),max_(0),err_(0){
    bounds_= usebounds_ = error_ = false;
    if (check_bounds(min, max)){
      min_ = min;
      max_ = max;
      bounds_ = true;
    }
    make_str();
  }

  //! Standard constructor with value, bounds and error
  /*!
   * Standard constructor with value, bounds and error provided. Creates
   * parameter with the given information. If a check for valid bounds
   * fails, just value and error are used.
   * \param value input value of the parameter
   * \param min input lower bound
   * \param max input upper bound
   * \param error input error of the parameter
   * \sa make_str(), check_bounds()
  */
  PWAGenericPar(const T value, const T min, const T max, const T error)
  :val_(value),min_(0),max_(0),err_(error){
    bounds_= usebounds_ = false;
    error_ = true;
    if (check_bounds(min, max)){
      min_ = min;
      max_ = max;
      bounds_ = true;
    }
    make_str();
  }

  //! Copy constructor using = operator
  /*!
   * Simple copy constructor using the = operator. As this operator is not
   * overloaded in this class, c++ will copy every member variable. As this
   * is a container class, this should be fine.
   * \param in input PWAParameter which variables will be copied
  */
  PWAGenericPar(const PWAGenericPar<T>& in){
    *this = in;
  }

  //! Empty Destructor
  /*!
   * There is nothing to destroy :(
  */
  virtual ~PWAGenericPar() { /* nothing */	}

  //! Getter for value of parameter cast to double
  virtual const inline double GetValue() const {return (double)val_;}
  //! Getter for lower bound of parameter cast to double
  virtual const inline double GetMinValue() const {return (double)min_;}
  //! Getter for upper bound of parameter cast to double
  virtual const inline double GetMaxValue() const {return (double)max_;}
  //! Getter for error of parameter cast to double
  virtual const inline double GetError() const {return (double)err_;}

  //! Getter for value of parameter
  virtual const inline T GetTValue() const {return val_;}
  //! Getter for lower bound of parameter
  virtual const inline T GetTMinValue() const {return min_;}
  //! Getter for upper bound of parameter
  virtual const inline T GetTMaxValue() const {return max_;}
  //! Getter for error of parameter
  virtual const inline T GetTError() const {return err_;}

  //! Setter for value of parameter cast from double
  virtual void SetValue(const double inVal) {val_ = (double)inVal; make_str();}
  //! Setter for error of parameter cast from double
  virtual void SetError(const double inErr) {err_ = (double)inErr; error_ = true; make_str();}
  //! Setter for bounds of parameter cast to double
  virtual const bool SetMinMax(const double inMin, const double inMax){
    bool valid = check_bounds((T)inMin, (T)inMax);
    if(valid){
        min_ = (T)inMin;
        max_ = (T)inMax;
        bounds_ = true;
        make_str();
    }
    return valid;
  }

  //! Setter for value of parameter
  virtual void SetTValue(const T value) {val_ = value; make_str();}
  //! Setter for error of parameter
  virtual void SetTError(const T error) {err_ = error; error_ = true; make_str();}

  //! Setter for bounds
  /*!
   * Setter for bounds of the parameter. If a check for valid bounds fails, it
   * returns false and nothing changes. This means if the bounds are invalid
   * the parameter maintains its old bounds if it had some.
   * \param min input lower bound
   * \param max input upper bound
   * \return bool if successful (re)set the bounds
   * \sa check_bounds()
  */
  virtual const bool SetTMinMax(const T min, const T max) {
    bool valid = check_bounds(min, max);
    if(valid){
        min_ = min;
        max_ = max;
        bounds_ = true;
        make_str();
    }
    return valid;
  }

  //! Setter for lower bound
  /*!
   * Setter for lower bound of the parameter. If a check for valid bounds
   * fails, it returns false and nothing changes. This means if the lower
   * bound is invalid the parameter maintains its old bounds if it had some.
   * \param min input lower bound
   * \return bool if successful (re)set lower bound
   * \sa check_bounds()
  */
  virtual const bool SetTMinValue(const T min) {
    bool valid = check_bounds(min, max_);
    if(valid){
        min_ = min;
        bounds_ = true;
        make_str();
    }
    return valid;
  }

  //! Setter for upper bound
  /*!
   * Setter for upper bound of the parameter. If a check for valid bounds
   * fails, it returns false and nothing changes. This means if the upper
   * bound is invalid the parameter maintains its old bounds if it had some.
   * \param max input upper bound
   * \return bool if successful (re)set upper bound
   * \sa check_bounds()
  */
  virtual const bool SetTMaxValue(const T max) {
    bool valid = check_bounds(min_, max);
    if(valid){
        max_ = max;
        bounds_ = true;
        make_str();
    }
    return valid;
  }

  //! A public function returning a string naming its type
  /*!
   * This function is used to get the type of the implementation of this
   * general parameter interface. Beware: using typid gives compiler-dependent
   * output!
   * \sa operator<<, to_str(), make_str()
  */
  virtual const std::string type(){
    return (std::string) typeid(T).name();
  }

protected:
  T val_, min_, max_, err_; /*!< Containers of parameter information */

  //! A protected function to check if bounds are valid
  /*!
   * This function checks if the bounds of the parameter are valid:
   * Upper bound should be larger then lower bound and the value
   * should be inside of the bounds.
   * \param max upper bound to check
   * \param min lower bound to check
   * \return bool if bounds are valid
   * \sa PWAParameter(const T value, const T min, const T max)
   * \sa PWAParameter(const T value, const T min, const T max, const T error)
   * \sa SetMinMax(), SetMinValue(), SetMaxValue()
  */
  bool check_bounds(const T min, const T max){
    if( (max > min) && (max >= val_) && (min <= val_))
      return true;
    return false;
  }

  //! A protected function which creates an output string for printing
  /*!
   * This function uses all available information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str(), type()
  */
  virtual void make_str() {
    std::stringstream oss;
    oss << "\t Val = " << val_;
    if(bounds_)
      oss << "\t  Min-Max = " << min_ << " to " << max_;
    if(error_)
      oss << "\t  Err = " << err_;
    oss << "\t Type = " << type();
    out_ = oss.str();
  }

};

#endif
