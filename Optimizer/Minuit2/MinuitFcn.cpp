#include "Optimizer/Minuit2/MinuitFcn.hpp"
#include "Optimizer/ControlParameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"
//#include "ErrLogger/ErrLogger.hh"
#include <cassert>
#include <memory>
#include <iostream>

using namespace ROOT::Minuit2;

MinuitFcn::MinuitFcn(std::shared_ptr<ControlParameter> myData) :
  _myDataPtr(myData){
  if (0==_myDataPtr) {
    //Alert << "Data pointer is 0 !!!!" << endmsg;
      std::cout << "Data pointer is 0 !!!!" << std::endl; //TODO exception
    exit(1);
  }
}

MinuitFcn::~MinuitFcn(){
  //std::cout << "~MinuitFcn: I'll be back" << std::endl;
}

double MinuitFcn::operator()(const std::vector<double>& x) const{
  ParameterList par;
  for(unsigned int i=0; i<x.size(); i++){
    std::cout << x[i] << " ";
    par.AddParameter(DoubleParameter(x[i]));
  }
  double result=_myDataPtr->controlParameter(par);
  std::cout << std::endl << "current minimized value:\t"<< result << std::endl;

  return result;
}

double MinuitFcn::Up() const{
return 0.5; //TODO: Setter, LH 0.5, Chi2 1.
}



