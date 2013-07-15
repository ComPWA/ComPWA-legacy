//! Wrapper of the Minuit2 Optimizer library.
/*! \class MinuitIF
 * @file MinuitIF.hpp
 * This class provides a wrapper around the Minuit2 library. It fulfills the
 * Optimizer interface to be easily adapted to other modules. The data needs to
 * be provided with the ControlParameter interface.
*/

#ifndef _OIFMINUIT_HPP
#define _OIFMINUIT_HPP

#include <vector>
//#include <boost/shared_ptr.hpp>
#include <memory>

#include "Optimizer/ControlParameter.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/ParameterList.hpp"
#include "Optimizer/Minuit2/MinuitFcn.hpp"

using namespace ROOT::Minuit2;

class MinuitIF : public Optimizer {

public:
  /// Default Constructor (0x0)
  MinuitIF(std::shared_ptr<ControlParameter> theData);
  virtual const double exec(ParameterList& par);

  /** Destructor */
  virtual ~MinuitIF();

 protected:

 private:
  MinuitFcn _myFcn;
 // vector<string> paramNames;
};

#endif /* _OIFMinuit_HPP */
