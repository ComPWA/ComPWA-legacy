#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include "PIFBW.hpp"
#include "PWAParameter.hpp"
#include "PWAGenericPar.hpp"

using namespace std;

PIFBW::PIFBW() {
  //_myFcn = new MIMinuitFcn(theData);
}

PIFBW::~PIFBW()
{
  //delete _myFcn;
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
  return true;
}

const double PIFBW::BreitWigner(double x, double M, double T){
  double denom=(x*x-M*M)*(x*x-M*M)+M*M*T*T;

  return 1./denom;
}
