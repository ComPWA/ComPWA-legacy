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

  PWAParameter() {
    val_ = 0;
    min_ = 0;
    max_ = 0;
    err_ = 0;
    make_str();
  }

  PWAParameter(const T value, const T min, const T max, const T error){
    val_ = value;
    min_ = min;
    max_ = max;
    err_ = error;
    make_str();
  }

  PWAParameter(const PWAParameter<T>& in){
    *this = in;
   /* val_ = in.GetValue();
    min_ = in.GetMinValue();
    max_ = in.GetMaxValue();
    err_ = in.GetError();*/
  }

  virtual ~PWAParameter() { /* nothing */	}

  virtual const inline T GetValue() const {return val_;}
  virtual const inline T GetMinValue() const {return min_;}
  virtual const inline T GetMaxValue() const {return max_;}
  virtual const inline T GetError() const {return err_;}

  virtual inline void SetValue(const T value) {val_ = value; make_str();}
  virtual inline void SetMinValue(const T min) {min_ = min; make_str();}
  virtual inline void SetMaxValue(const T max) {max_ = max; make_str();}
  virtual inline void SetError(const T error) {err_ = error; make_str();}

  std::string const& to_str() const{
    return out;
  }

protected:
  T val_, min_, max_, err_;
  std::string out;
  //TODO: other parameter info?

  void make_str() {
    std::stringstream oss;
    oss << "\t Val = " << val_ << "\t  Min-Max = " << min_ << " to " << max_ << "\t  Err = " << err_;
    out = oss.str();
  }

  template <class U>
  friend std::ostream & operator<<(std::ostream &os, const PWAParameter<U>& p);

};

template <class U>
std::ostream & operator<<(std::ostream &os, const PWAParameter<U>& p){
  return os << p.to_str();
}

#endif
