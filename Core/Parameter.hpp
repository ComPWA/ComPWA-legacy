//! Base class for internal parameter.
/*! \class DoubleParameter
 * \class IntegerParameter
 * \class BoolParameter
 * @file Parameter.hpp
 * This class defines the internal container of a fit parameter.
 * A parameter consists of a value with optional bounds and error.
*/

#ifndef _PARAMETER_HPP_
#define _PARAMETER_HPP_

#include <iostream>
#include <string>
#include <sstream>

#include "Core/Exceptions.hpp"

class DoubleParameter
{

public:

  //! Standard constructor without information
   /*!
    * Standard constructor with no information provided. Creates parameter
    * with value 0 but without bounds or an error.
    * \sa make_str()
   */
  DoubleParameter():val_(0),min_(0),max_(0),err_(0) {
     bounds_= usebounds_ = error_ = fixed_ = false;
     make_str();
   }

   //! Standard constructor with a value
   /*!
    * Standard constructor with just a value provided. Creates parameter
    * with given value but without bounds or an error.
    * \param value input value of the parameter
    * \sa make_str()
   */
  DoubleParameter(const double value):val_(value),min_(0),max_(0),err_(0){
     bounds_= usebounds_ = error_ = fixed_ = false;
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
  DoubleParameter(const double value, const double error)
   :val_(value),min_(0),max_(0),err_(error){
     bounds_= usebounds_ = fixed_ = false;
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
  DoubleParameter(const double value, const double min, const double max)
   :val_(value),min_(0),max_(0),err_(0){
     bounds_= usebounds_ = error_ = fixed_ = false;
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
  DoubleParameter(const double value, const double min, const double max, const double error)
   :val_(value),min_(0),max_(0),err_(error){
     bounds_= usebounds_ = fixed_ = false;
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
  DoubleParameter(const DoubleParameter& in){
     *this = in;
   }

  //! Empty Destructor
  /*!
   * There is nothing to destroy :(
  */
  virtual ~DoubleParameter() { /* nothing */	}

  //! Check if parameter has bounds
  virtual const inline bool HasBounds() const {return bounds_;}
  //! Check if bounds should be used
  virtual const inline bool UseBounds() const {if(bounds_)return usebounds_; return false;}
  //! Check if parameter has an error
  virtual const inline bool HasError() const {return error_;}
  //! Check if parameter is fixed
  virtual const inline bool IsFixed() const {return fixed_;}

  //! Getter for value of parameter
  virtual const inline double GetValue() const {return val_;}
  //! Getter for lower bound of parameter
  virtual const inline double GetMinValue() const {return min_;}
  //! Getter for upper bound of parameter
  virtual const inline double GetMaxValue() const {return max_;}
  //! Getter for error of parameter
  virtual const inline double GetError() const {return err_;}

  //! Setter for value of parameter
  virtual void SetValue(const double inVal) {
    if(fixed_){
        throw ParameterFixed();
        return;
    }
    val_ = inVal;
    make_str();
  }
  //! Setter for error of parameter
  virtual void SetError(const double inErr) {err_ = inErr; error_ = true; make_str();}
  //! Setter for bounds of parameter
  virtual const bool SetMinMax(const double inMin, const double inMax){
    bool valid = check_bounds(inMin, inMax);
    if(valid){
        min_ = inMin;
        max_ = inMax;
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
  virtual const bool SetMinValue(const double min) {
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
  virtual const bool SetMaxValue(const double max) {
    bool valid = check_bounds(min_, max);
    if(valid){
        max_ = max;
        bounds_ = true;
        make_str();
    }
    return valid;
  }

  //! Set if bounds should be used
  virtual const inline void UseBounds(const bool use) {usebounds_=use;}
  //! Call to fix parameter
  virtual const inline void SetParameterFixed() {fixed_=true;}
  //! Call to free parameter
  virtual const inline void SetParameterFree() {fixed_=false;}
  //! Set parameter free or fixed
  virtual const inline void FixParameter(const bool fixed) {fixed_=fixed;}

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
   * \sa operator<<, to_str(), make_str()
  */
  virtual const std::string type(){
    return "double";
  }

protected:
  std::string out_; /*!< Output string to print information */
  bool bounds_; /*!< Are valid bounds defined for this parameter? */
  bool error_; /*!< Is an error defined for this parameter? */
  bool usebounds_; /*!< Do you want to restrict your parameter? */
  bool fixed_; /*!< Do you want to keep parameter fixed? */
  double val_, min_, max_, err_; /*!< Containers of parameter information */

  //! A protected function to check if bounds are valid
  /*!
   * This function checks if the bounds of the parameter are valid:
   * Upper bound should be larger then lower bound and the value
   * should be inside of the bounds.
   * \param max upper bound to check
   * \param min lower bound to check
   * \return bool if bounds are valid
   * \sa Parameter(const double value, const double min, const double max)
   * \sa Parameter(const double value, const double min, const double max, const double error)
   * \sa SetMinMax(), SetMinValue(), SetMaxValue()
  */
  bool check_bounds(const double min, const double max){
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

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type. The definition
   * of this class has to be outside the namespace of the class.
   * \sa make_str(), to_str()
  */
  friend std::ostream & operator<<(std::ostream &os, const DoubleParameter& p);

};

inline std::ostream & operator<<(std::ostream &os, const DoubleParameter& p){
  return os << p.to_str();
}

class IntegerParameter
{

public:

  //! Standard constructor without information
   /*!
    * Standard constructor with no information provided. Creates parameter
    * with value 0 but without bounds or an error.
    * \sa make_str()
   */
  IntegerParameter():val_(0),min_(0),max_(0),err_(0) {
     bounds_= usebounds_ = error_ = fixed_ = false;
     make_str();
   }

   //! Standard constructor with a value
   /*!
    * Standard constructor with just a value provided. Creates parameter
    * with given value but without bounds or an error.
    * \param value input value of the parameter
    * \sa make_str()
   */
  IntegerParameter(const int value):val_(value),min_(0),max_(0),err_(0){
     bounds_= usebounds_ = error_ = fixed_ = false;
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
  IntegerParameter(const int value, const int error)
   :val_(value),min_(0),max_(0),err_(error){
     bounds_= usebounds_ = fixed_ = false;
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
  IntegerParameter(const int value, const int min, const int max)
   :val_(value),min_(0),max_(0),err_(0){
     bounds_= usebounds_ = error_ = fixed_ = false;
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
  IntegerParameter(const int value, const int min, const int max, const int error)
   :val_(value),min_(0),max_(0),err_(error){
     bounds_= usebounds_ = fixed_ = false;
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
  IntegerParameter(const IntegerParameter& in){
     *this = in;
   }

  //! Empty Destructor
  /*!
   * There is nothing to destroy :(
  */
  virtual ~IntegerParameter() { /* nothing */    }

  //! Check if parameter has bounds
  virtual const inline bool HasBounds() const {return bounds_;}
  //! Check if bounds should be used
  virtual const inline bool UseBounds() const {if(bounds_)return usebounds_; return false;}
  //! Check if parameter has an error
  virtual const inline bool HasError() const {return error_;}
  //! Check if parameter is fixed
  virtual const inline bool IsFixed() const {return fixed_;}

  //! Getter for value of parameter
  virtual const inline int GetValue() const {return val_;}
  //! Getter for lower bound of parameter
  virtual const inline int GetMinValue() const {return min_;}
  //! Getter for upper bound of parameter
  virtual const inline int GetMaxValue() const {return max_;}
  //! Getter for error of parameter
  virtual const inline int GetError() const {return err_;}

  //! Setter for value of parameter
  virtual void SetValue(const int inVal) {
    if(fixed_){
        throw ParameterFixed();
        return;
    }
    val_ = inVal;
    make_str();
  }
  //! Setter for error of parameter
  virtual void SetError(const int inErr) {err_ = inErr; error_ = true; make_str();}
  //! Setter for bounds of parameter
  virtual const bool SetMinMax(const int inMin, const int inMax){
    bool valid = check_bounds(inMin, inMax);
    if(valid){
        min_ = inMin;
        max_ = inMax;
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
  virtual const bool SetMinValue(const int min) {
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
  virtual const bool SetMaxValue(const int max) {
    bool valid = check_bounds(min_, max);
    if(valid){
        max_ = max;
        bounds_ = true;
        make_str();
    }
    return valid;
  }

  //! Set if bounds should be used
  virtual const inline void UseBounds(const bool use) {usebounds_=use;}
  //! Call to fix parameter
  virtual const inline void SetParameterFixed() {fixed_=true;}
  //! Call to free parameter
  virtual const inline void SetParameterFree() {fixed_=false;}
  //! Set parameter free or fixed
  virtual const inline void FixParameter(const bool fixed) {fixed_=fixed;}

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
   * \sa operator<<, to_str(), make_str()
  */
  virtual const std::string type(){
    return "integer";
  }

protected:
  std::string out_; /*!< Output string to print information */
  bool bounds_; /*!< Are valid bounds defined for this parameter? */
  bool error_; /*!< Is an error defined for this parameter? */
  bool usebounds_; /*!< Do you want to restrict your parameter? */
  bool fixed_; /*!< Do you want to keep parameter fixed? */
  int val_, min_, max_, err_; /*!< Containers of parameter information */

  //! A protected function to check if bounds are valid
  /*!
   * This function checks if the bounds of the parameter are valid:
   * Upper bound should be larger then lower bound and the value
   * should be inside of the bounds.
   * \param max upper bound to check
   * \param min lower bound to check
   * \return bool if bounds are valid
   * \sa Parameter(const double value, const double min, const double max)
   * \sa Parameter(const double value, const double min, const double max, const double error)
   * \sa SetMinMax(), SetMinValue(), SetMaxValue()
  */
  bool check_bounds(const int min, const int max){
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

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type. The definition
   * of this class has to be outside the namespace of the class.
   * \sa make_str(), to_str()
  */
  friend std::ostream & operator<<(std::ostream &os, const IntegerParameter& p);

};

inline std::ostream & operator<<(std::ostream &os, const IntegerParameter& p){
  return os << p.to_str();
}

class BoolParameter
{

public:

  //! Standard constructor without information
   /*!
    * Standard constructor with no information provided. Creates parameter
    * with value 0 but without bounds or an error.
    * \sa make_str()
   */
  BoolParameter():val_(0),err_(0) {
     error_ = fixed_ = false;
     make_str();
   }

   //! Standard constructor with a value
   /*!
    * Standard constructor with just a value provided. Creates parameter
    * with given value but without bounds or an error.
    * \param value input value of the parameter
    * \sa make_str()
   */
  BoolParameter(const bool value):val_(value),err_(0){
     error_ = fixed_ = false;
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
  BoolParameter(const bool value, const bool error)
   :val_(value),err_(error){
     fixed_ = false;
     error_ = true;
     make_str();
   }

   //! Copy constructor using = operator
   /*!
    * Simple copy constructor using the = operator. As this operator is not
    * overloaded in this class, c++ will copy every member variable. As this
    * is a container class, this should be fine.
    * \param in input PWAParameter which variables will be copied
   */
  BoolParameter(const BoolParameter& in){
     *this = in;
   }

  //! Empty Destructor
  /*!
   * There is nothing to destroy :(
  */
  virtual ~BoolParameter() { /* nothing */    }

  //! Check if parameter has an error
  virtual const inline bool HasError() const {return error_;}
  //! Check if parameter is fixed
  virtual const inline bool IsFixed() const {return fixed_;}

  //! Getter for value of parameter
  virtual const inline bool GetValue() const {return val_;}
  //! Getter for error of parameter
  virtual const inline bool GetError() const {return err_;}

  //! Setter for value of parameter
  virtual void SetValue(const bool inVal) {
    if(fixed_){
        throw ParameterFixed();
        return;
    }
    val_ = inVal;
    make_str();
  }
  //! Setter for error of parameter
  virtual void SetError(const bool inErr) {err_ = inErr; error_ = true; make_str();}

  //! Call to fix parameter
  virtual const inline void SetParameterFixed() {fixed_=true;}
  //! Call to free parameter
  virtual const inline void SetParameterFree() {fixed_=false;}
  //! Set parameter free or fixed
  virtual const inline void FixParameter(const bool fixed) {fixed_=fixed;}


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
   * \sa operator<<, to_str(), make_str()
  */
  virtual const std::string type(){
    return "boolean";
  }

protected:
  std::string out_; /*!< Output string to print information */
  bool error_; /*!< Is an error defined for this parameter? */
  bool usebounds_; /*!< Do you want to restrict your parameter? */
  bool fixed_; /*!< Do you want to keep parameter fixed? */
  int val_, err_; /*!< Containers of parameter information */

  //! A protected function which creates an output string for printing
  /*!
   * This function uses all available information about the parameter
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str(), type()
  */
  virtual void make_str() {
    std::stringstream oss;
    oss << "\t Val = " << val_;
    if(error_)
      oss << "\t  Err = " << err_;
    oss << "\t Type = " << type();
    out_ = oss.str();
  }

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type. The definition
   * of this class has to be outside the namespace of the class.
   * \sa make_str(), to_str()
  */
  friend std::ostream & operator<<(std::ostream &os, const BoolParameter& p);

};

inline std::ostream & operator<<(std::ostream &os, const BoolParameter& p){
  return os << p.to_str();
}

#endif
