#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <memory>
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"
#include "OIFMinuit.hpp"
//#include "Examples/Tutorial/FitIF/MinimizerInterface/MIBase.hh"

using namespace ROOT::Minuit2;

OIFMinuit::OIFMinuit(std::shared_ptr<OIFData> theData) : _myFcn(theData){
  //_myFcn = new MIMinuitFcn(theData);
}

OIFMinuit::~OIFMinuit(){
  //delete _myFcn;
}

const double OIFMinuit::exec(unsigned int Npar,  double* par,  double* min, double* max, double* err){

  std::string s;
  std::stringstream out;
  
  MnUserParameters upar;
  for(unsigned int i=0; i<Npar; ++i){
    out.str("");
    out << i;
    s = out.str();
    upar.Add(s, par[i], err[i], max[i], min[i]);
  }

  MnMigrad migrad(_myFcn, upar);
 // Info <<"start migrad "<< endmsg;
  FunctionMinimum minMin = migrad();

 if(!minMin.IsValid()) {
   //try with higher strategy
 //  Info <<"FM is invalid, try with strategy = 2."<< endmsg;
   MnMigrad migrad2(_myFcn, minMin.UserState(), MnStrategy(2));
   minMin = migrad2();
 }
 
  //save minized values
  for(unsigned int i=0; i<Npar; ++i){
    out.str("");
    out << i;
    s = out.str();
    par[i]=minMin.UserState().Value(s); err[i]=minMin.UserState().Error(s);
  }

  return minMin.Fval();
}
