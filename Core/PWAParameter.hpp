//! Internal container representing a parameter.
/*! \class PWAParameter
 * @file PWAParameter.hpp
 * This class provides a internal container for information of a fit parameter.
 * The class provides the value, bounds and error of a parameter.
*/

#ifndef _PWAPARAMETER_HPP_
#define _PWAPARAMETER_HPP_

#include <iostream>
#include <string>
#include <sstream>

template <class T>
class PWAParameter
{

public:

  PWAParameter():val_(0),min_(0),max_(0),err_(0),bounds_(false),error_(false) {
    make_str();
  }

  PWAParameter(const T value):val_(value),min_(0),max_(0),err_(0),bounds_(false),error_(false){
    make_str();
  }

  PWAParameter(const T value, const T error)
  :val_(value),min_(0),max_(0),err_(error),bounds_(false),error_(true){
    make_str();
  }

  PWAParameter(const T value, const T min, const T max)
  :val_(value),min_(min),max_(max),err_(0),bounds_(true),error_(false){
    check_bounds();
    make_str();
  }

  PWAParameter(const T value, const T min, const T max, const T error)
  :val_(value),min_(min),max_(max),err_(error),bounds_(true),error_(true){
    check_bounds();
    make_str();
  }

  PWAParameter(const PWAParameter<T>& in){
    *this = in;
  }

  virtual ~PWAParameter() { /* nothing */	}

  virtual const inline bool HasBounds() const {return bounds_;}
  virtual const inline bool HasError() const {return error_;}

  virtual const inline T GetValue() const {return val_;}
  virtual const inline T GetMinValue() const {return min_;}
  virtual const inline T GetMaxValue() const {return max_;}
  virtual const inline T GetError() const {return err_;}

  virtual void SetValue(const T value) {val_ = value; make_str();}
  virtual void SetError(const T error) {err_ = error; make_str();}
  virtual const bool SetMinMax(const T min, const T max) {
    min_ = min; max_ = max;
    bool worked = check_bounds();
    make_str();
    return worked;
  }
  virtual const bool SetMinValue(const T min) {
    min_ = min;
    bool worked = check_bounds();
    make_str();
    return worked;
  }
  virtual const bool SetMaxValue(const T max) {
    max_ = max;
    bool worked = check_bounds();
    make_str();
    return worked;
  }

  std::string const& to_str() const{
    return out_;
  }

protected:
  T val_, min_, max_, err_;
  std::string out_;
  bool bounds_;
  bool error_;

  //TODO: other parameter info?

  bool check_bounds(){
    if( max_ <= min_){
      max_=min_ = 0;
      bounds_ = false;
      //TODO: exception, message, something
    } else {
      bounds_ = true;
    }
    return bounds_;
  }

  void make_str() {
    std::stringstream oss;
    oss << "\t Val = " << val_;
    if(bounds_)
      oss << "\t  Min-Max = " << min_ << " to " << max_;
    if(error_)
      oss << "\t  Err = " << err_;
    out_ = oss.str();
  }

  template <class U>
  friend std::ostream & operator<<(std::ostream &os, const PWAParameter<U>& p);

};

template <class U>
std::ostream & operator<<(std::ostream &os, const PWAParameter<U>& p){
  return os << p.to_str();
}

#endif
