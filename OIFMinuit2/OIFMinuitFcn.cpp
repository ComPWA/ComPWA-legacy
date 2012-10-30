#include "OIFMinuitFcn.hpp"
#include "OIFData.hpp"
//#include "ErrLogger/ErrLogger.hh"
#include <cassert>
#include <memory>
#include <iostream>

using namespace ROOT::Minuit2;
using namespace std;

OIFMinuitFcn::OIFMinuitFcn(shared_ptr<OIFData> myData) :
  _myDataPtr(myData)
{
  if (0==_myDataPtr) {
    //Alert << "Data pointer is 0 !!!!" << endmsg;
	cout << "Data pointer is 0 !!!!" << endl;
    exit(1);
  }
}

OIFMinuitFcn::~OIFMinuitFcn()
{
}

double OIFMinuitFcn::operator()(const vector<double>& par) const
{
  double result=_myDataPtr->controlParameter(par);
  //DebugMsg << "current minimized value:\t"<< result << endmsg;
  cout << "current minimized value:\t"<< result << endl;
  return result;
}

double OIFMinuitFcn::Up() const
{
return 1.;
}



