//! Minuit2 function to be optimized.
/*! \class OIFMinuitFcn
 * @file OIFMinuitFcn.hpp
 * Based on the Minuit2 FcnBase. This class uses the OIFData interface for the
 * optimization.
*/

#ifndef _OIFMinuitFcn_HPP
#define _OIFMinuitFcn_HPP

#include <vector>
#include <memory>
//#include <boost/shared_ptr.hpp>
//#include <cassert>
#include "Minuit2/FCNBase.h"
#include "OIFData.hpp"
#include "PWAParameter.hpp"

class OIFData;

namespace ROOT {

   namespace Minuit2 {
class OIFMinuitFcn : public FCNBase {

public:

  OIFMinuitFcn(std::shared_ptr<OIFData> theData);
  virtual ~OIFMinuitFcn();

  double operator()(const std::vector<double>& x) const;

  double Up() const;

private:
  std::shared_ptr<OIFData> _myDataPtr;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif 
