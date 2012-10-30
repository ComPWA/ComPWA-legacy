#ifndef _OIFMINUIT_HPP
#define _OIFMINUIT_HPP

#include <vector>
//#include <boost/shared_ptr.hpp>
#include <memory>
#include "OIFData.hpp"
#include "OIFBase.hpp"
#include "OIFMinuitFcn.hpp"

using namespace ROOT::Minuit2;
using namespace std;

class OIFMinuit : public OIFBase {

public:
  /// Default Constructor (0x0)
	OIFMinuit(shared_ptr<OIFData> theData);
  virtual const double exec(unsigned int Npar, double* par,  double* min, double* max, double* err); 

  /** Destructor */
  virtual ~OIFMinuit();

 protected:

 private:
  OIFMinuitFcn _myFcn;
 // vector<string> paramNames;
};

#endif /* _OIFMinuit_HPP */
