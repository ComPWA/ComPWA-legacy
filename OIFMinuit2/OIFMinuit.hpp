//! Wrapper of the Minuit2 Optimizer library.
/*! \class OIFMinuit
 * @file OIFMinuit.hpp
 * This class provides a wrapper around the Minuit2 library. It fulfills the
 * Optimizer interface to be easily adapted to other modules. The data needs to
 * be provided with the ControlParameter interface.
*/

#ifndef _OIFMINUIT_HPP
#define _OIFMINUIT_HPP

#include <vector>
//#include <boost/shared_ptr.hpp>
#include <memory>

#include "ControlParameter.hpp"
#include "Optimizer.hpp"
#include "PWAParameter.hpp"
#include "OIFMinuitFcn.hpp"

using namespace ROOT::Minuit2;

class OIFMinuit : public Optimizer {

public:
  /// Default Constructor (0x0)
  OIFMinuit(std::shared_ptr<ControlParameter> theData);
  virtual const double exec(std::vector<std::shared_ptr<PWAParameter> >& par);

  /** Destructor */
  virtual ~OIFMinuit();

 protected:

 private:
  OIFMinuitFcn _myFcn;
 // vector<string> paramNames;
};

#endif /* _OIFMinuit_HPP */
