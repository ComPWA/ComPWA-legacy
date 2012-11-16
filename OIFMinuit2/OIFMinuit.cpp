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
#include "PWAParameter.hpp"

using namespace ROOT::Minuit2;

OIFMinuit::OIFMinuit(std::shared_ptr<OIFData> theData) : _myFcn(theData){
  //_myFcn = new MIMinuitFcn(theData);
}

OIFMinuit::~OIFMinuit(){
  //std::cout << "OIFMinuit::~OIFMinuit: I'll be back" << std::endl;
  //delete _myFcn;
}

const double OIFMinuit::exec(std::vector<PWAParameter<double> >& par){
  std::string s;
  std::stringstream out;

  MnUserParameters upar;
  for(unsigned int i=0; i<par.size(); ++i){
    out.str("");
    out << i;
    s = out.str();

    //use as much information as possible: (just bounds but no error not supported by minuit)
    if( par[i].HasBounds() && par[i].HasError() )
      upar.Add(s, par[i].GetValue(), par[i].GetError(), par[i].GetMaxValue(), par[i].GetMinValue());
    else if( par[i].HasError() )
      upar.Add(s, par[i].GetValue(), par[i].GetError());
    else
      upar.Add(s, par[i].GetValue());
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
  for(unsigned int i=0; i<par.size(); ++i){
    out.str("");
    out << i;
    s = out.str();
    par[i].SetValue(minMin.UserState().Value(s));
    par[i].SetError(minMin.UserState().Error(s));
  }

  return minMin.Fval();
}

