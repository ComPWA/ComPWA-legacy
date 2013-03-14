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
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"

using namespace ROOT::Minuit2;

MinuitIF::MinuitIF(std::shared_ptr<ControlParameter> theData) : _myFcn(theData){
  //_myFcn = new MIMinuitFcn(theData);
}

MinuitIF::~MinuitIF(){
  //std::cout << "MinuitIF::~MinuitIF: I'll be back" << std::endl;
  //delete _myFcn;
}

const double MinuitIF::exec(ParameterList& par){
  std::string s;
  std::stringstream out;

  MnUserParameters upar;
  for(unsigned int i=0; i<par.GetNDouble(); ++i){ //only doubles for minuit
    out.str("");
    out << i;
    s = out.str();

    //use as much information as possible: (just bounds but no error not supported by minuit)
    //try{
      DoubleParameter& actPat = par.GetDoubleParameter(i);
    //}
    if( actPat.UseBounds() && actPat.HasError() )
      upar.Add(s, actPat.GetValue(), actPat.GetError(), actPat.GetMaxValue(), actPat.GetMinValue());
    else if( actPat.HasError() )
      upar.Add(s, actPat.GetValue(), actPat.GetError());
    else
      upar.Add(s, actPat.GetValue());
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

  //save minimzed values
  for(unsigned int i=0; i<par.GetNDouble(); ++i){
    out.str("");
    out << i;
    s = out.str();
    DoubleParameter& actPat = par.GetDoubleParameter(i);
    actPat.SetValue(minMin.UserState().Value(s));
    actPat.SetError(minMin.UserState().Error(s));
  }

  return minMin.Fval();
}

