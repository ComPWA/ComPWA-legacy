#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include "Physics/BreitWigner/BreitWigner.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"

using namespace std;

BreitWigner::BreitWigner(const double min, const double max):min_(min),max_(max) {
  //_myFcn = new MIMinuitFcn(theData);
}

BreitWigner::~BreitWigner()
{
  //delete _myFcn;
}

const double BreitWigner::integral(ParameterList& par){
  double integral = 0;
  unsigned int nSteps = 1000000;
  double step = (max_-min_)/(double)nSteps;

  //TODO: try exceptions, parameter find via name?
  integral += step*BreitWignerValue(min_, par.GetDoubleParameter(0)->GetValue(), par.GetDoubleParameter(1)->GetValue())/2.;
  for(unsigned int k=1; k<nSteps; k++){
    integral += step*BreitWignerValue((min_+k*step), par.GetDoubleParameter(0)->GetValue(), par.GetDoubleParameter(1)->GetValue());
  }
  integral += step*BreitWignerValue(max_, par.GetDoubleParameter(0)->GetValue(), par.GetDoubleParameter(1)->GetValue())/2.;

  return integral;
}

const double BreitWigner::volume(){
  return 0;
}

const double BreitWigner::drawInt(double* x, double *p){
  return p[2]*BreitWignerValue(x[0], p[0], p[1]);
}

const ParameterList BreitWigner::intensity(double x, double M, double T){
  ParameterList result;
  result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BreitWignerResult",BreitWignerValue(x, M, T))));
  return result;
}

const ParameterList BreitWigner::intensity(std::vector<double>& x, ParameterList& par){
  ParameterList result;
  double val = BreitWignerValue(x.at(0), par.GetDoubleParameter(0)->GetValue(), par.GetDoubleParameter(1)->GetValue());
  result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BreitWignerResult",val)));
  return result;
}

const bool BreitWigner::fillStartParVec(ParameterList& outPar){
  if(outPar.GetNParameter())
    return false; //already filled ,TODO: exception?

  outPar.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BreitWignerPosition",1.5, 0.5, 2.5, 0.1)));
  outPar.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BreitWignerWidth",0.3, 0.1, 0.5, 0.05)));

  return true;
}

const double BreitWigner::BreitWignerValue(double x, double M, double T){
  double denom=(x*x-M*M)*(x*x-M*M)+M*M*T*T;

  return 1./denom;
}
