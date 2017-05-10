/* File : CoreToScript.i */
%module core
%{
/* Put headers and other declarations here */
#include "Parameter.hpp"
#include "Particle.hpp"
%}

%include "Particle.hpp"

namespace ComPWA {
class DoubleParameter : public AbsParameter {

public:
  DoubleParameter(std::string inName = "");

  DoubleParameter(std::string inName, const double value);

  DoubleParameter(std::string inName, const double value, const double error);

  DoubleParameter(std::string inName, const double value, const double min,
                  const double max);

  DoubleParameter(std::string inName, const double value, const double min,
                  const double max, const double error);
                  
  DoubleParameter(const DoubleParameter &in);

  virtual ~DoubleParameter();

  virtual inline bool HasBounds();

  virtual inline bool UseBounds();

  virtual inline bool IsFixed();

  virtual inline void SetUseBounds(const bool use);

  virtual inline void SetParameterFixed() ;

  virtual inline void SetParameterFree() ;

  virtual inline void FixParameter(const bool fixed) ;
  
  virtual void UpdateParameter(std::shared_ptr<DoubleParameter> newPar);

  virtual inline double GetValue();

  virtual inline double GetRoundedValue();

  virtual inline double GetMinValue();

  virtual inline double GetMaxValue();

  virtual std::complex<double> getNodeValue();

  virtual void SetValue(const double inVal);

  virtual void SetMinMax(const double min, const double max);

  virtual void SetMinValue(const double min);

  virtual void SetMaxValue(const double max);

  virtual inline bool HasError();

  virtual ErrorType GetErrorType();

  virtual double GetError();

  virtual double GetRoundedError();

  virtual double GetErrorHigh();

  virtual double GetErrorLow();

  virtual void SetError(double errLow, double errHigh);

  virtual void SetError(double err) ;

};
}