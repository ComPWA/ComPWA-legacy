//! Internal container representing a parameter.
/*! \class PWAParameter
 * @file PWAParameter.hpp
 * This class provides a internal container for information of a fit parameter.
 * The class provides the value, bounds and error of a parameter.
*/

#ifndef _PWAPARAMETER_HPP_
#define _PWAPARAMETER_HPP_

template <class T>
class PWAParameter
{

public:

  PWAParameter() {
    val_ = 0;
    min_ = 0;
    max_ = 0;
    err_ = 0;
  }

  PWAParameter(const T value, const T min, const T max, const T error){
    val_ = value;
    min_ = min;
    max_ = max;
    err_ = error;
  }

  PWAParameter(const PWAParameter& in){
    val_ = in.GetValue();
    min_ = in.GetMinValue();
    max_ = in.GetMaxValue();
    err_ = in.GetError();
  }

  virtual ~PWAParameter() { /* nothing */	}

  virtual const inline T GetValue() const {return val_;}
  virtual const inline T GetMinValue() const {return min_;}
  virtual const inline T GetMaxValue() const {return max_;}
  virtual const inline T GetError() const {return err_;}

  virtual const inline void SetValue(const T value) const {val_ = value;}
  virtual const inline void SetMinValue(const T min) const {min_ = min;}
  virtual const inline void SetMaxValue(const T max) const {max_ = max;}
  virtual const inline void SetError(const T error) const {err_ = error;}

protected:
  T val_, min_, max_, err_;
  //TODO: other parameter info?

};

#endif
