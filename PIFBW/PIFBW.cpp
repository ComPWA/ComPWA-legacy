#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include "PIFBW.hpp"
#include "PWAParameter.hpp"
#include "PWAGenericPar.hpp"

using namespace std;

PIFBW::PIFBW(const double min, const double max) : PIFBase(min, max) {
  //_myFcn = new MIMinuitFcn(theData);
}

PIFBW::~PIFBW()
{
  //delete _myFcn;
}

const double PIFBW::integral(std::vector<std::shared_ptr<PWAParameter> >& par){
  double integral = 0;
  unsigned int nSteps = 1000000;
  double step = (max_-min_)/(double)nSteps;

  integral += step*BreitWigner(min_, par.at(0)->GetValue(), par.at(1)->GetValue())/2.;
  for(unsigned int k=1; k<nSteps; k++){
    integral += step*BreitWigner((min_+k*step), par.at(0)->GetValue(), par.at(1)->GetValue());
  }
  integral += step*BreitWigner(max_, par.at(0)->GetValue(), par.at(1)->GetValue())/2.;

  return integral;
}
const double PIFBW::drawInt(double* x, double *p){
  return p[2]*BreitWigner(x[0], p[0], p[1]);
}

const double PIFBW::intensity(double x, double M, double T){
  return BreitWigner(x, M, T);
}

const double PIFBW::intensity(std::vector<double> x, std::vector<std::shared_ptr<PWAParameter> >& par){
  return BreitWigner(x.at(0), par.at(0)->GetValue(), par.at(1)->GetValue());
}

const bool PIFBW::fillStartParVec(std::vector<std::shared_ptr<PWAParameter> >& outPar){
  if(outPar.size())
    return false; //already filled
  outPar.push_back(shared_ptr<PWAParameter>(new PWAGenericPar<double>(1.5, 0.5, 2.5, 0.1)));
  outPar.push_back(shared_ptr<PWAParameter>(new PWAGenericPar<double>(0.3, 0.1, 0.5, 0.05)));
  //outPar.push_back(shared_ptr<PWAParameter>(new PWAGenericPar<double>(1.)));
  return true;
}

const double PIFBW::BreitWigner(double x, double M, double T){
  double denom=(x*x-M*M)*(x*x-M*M)+M*M*T*T;

  return 1./denom;
}
