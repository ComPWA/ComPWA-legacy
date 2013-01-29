#include "Optimizer/Minuit2/MinuitFcn.hpp"
#include "Core/ControlParameter.hpp"
#include "Core/PWAParameter.hpp"
#include "Core/PWAGenericPar.hpp"
//#include "ErrLogger/ErrLogger.hh"
#include <cassert>
#include <memory>
#include <iostream>

using namespace ROOT::Minuit2;

MinuitFcn::MinuitFcn(std::shared_ptr<ControlParameter> myData) :
  _myDataPtr(myData){
  if (0==_myDataPtr) {
    //Alert << "Data pointer is 0 !!!!" << endmsg;
      std::cout << "Data pointer is 0 !!!!" << std::endl;
    exit(1);
  }
}

MinuitFcn::~MinuitFcn(){
  //std::cout << "~MinuitFcn: I'll be back" << std::endl;
}

double MinuitFcn::operator()(const std::vector<double>& x) const{
  std::vector<std::shared_ptr<PWAParameter> > par;
  for(unsigned int i=0; i<x.size(); i++)
    par.push_back(std::shared_ptr<PWAGenericPar<double> >(new PWAGenericPar<double>(x[i])));
  double result=_myDataPtr->controlParameter(par);
  //DebugMsg << "current minimized value:\t"<< result << endmsg;
  std::cout << "current minimized value:\t"<< result << std::endl;
  return result;
}

double MinuitFcn::Up() const{
return 1.;
}



