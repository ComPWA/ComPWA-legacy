
#ifndef AMPFCN_HPP
#define AMPFCN_HPP

//#include <bits/shared_ptr.h>

#include "../Core/Amplitude.hpp"
#include "Estimator.hpp"

namespace ComPWA {

class ParameterList;

namespace Estimator {

class AmpFcn : public Estimator {

public:
  //	AmpFcn (std::shared_ptr<Amplitude> amp):_amp(amp){};
  AmpFcn(Amplitude *amp) : _amp(amp) { _amp->FillParameterList(_par); };
  //	virtual double controlParameter(std::vector<double> x, ParameterList&
  //minPar){
  //		return _amp->intensity(x,minPar);
  //	};
  virtual double controlParameter(ParameterList &x) {
    if (x.GetNParameter() != 2)
      return -999;
    std::vector<double> xx; // convert parameterList to std::vector
    xx.push_back(x.GetDoubleParameter(0)->GetValue());
    xx.push_back(x.GetDoubleParameter(1)->GetValue());
    //		xx.push_back(x.GetDoubleParameter(2).GetValue());
    //		for(unsigned int i=0;i< _par.GetNParameter();i++)
    //std::cout<<_par.GetParameterValue(i)<<std::endl;
    //		std::cout<<"===="<<std::endl;
    //		for(unsigned int i=0;i< xx.size();i++)
    //std::cout<<xx[i]<<std::endl;
    ParameterList result = _amp->intensity(xx);
    return ((-1) * result.GetDoubleParameter(0)->GetValue());
  };

protected:
private:
  //	std::shared_ptr<Amplitude> _amp;
  Amplitude *_amp;
  ParameterList _par;
};

} /* namespace Estimator */
} /* namespace ComPWA */

#endif
