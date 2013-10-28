//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Minuit2 function to be optimized.
/*! \class MinuitFcn
 * @file MinuitFcn.hpp
 * Based on the Minuit2 FcnBase. This class uses the ControlParameter interface for the
 * optimization.
*/

#ifndef _OIFMinuitFcn_HPP
#define _OIFMinuitFcn_HPP

#include <vector>
#include <memory>
//#include <boost/shared_ptr.hpp>
//#include <cassert>
#include "Minuit2/FCNBase.h"
#include "Optimizer/ControlParameter.hpp"
#include "Core/ParameterList.hpp"

class ControlParameter;

namespace ROOT {

   namespace Minuit2 {
class MinuitFcn : public FCNBase {

public:

  MinuitFcn(std::shared_ptr<ControlParameter> theData, ParameterList& parList);
  virtual ~MinuitFcn();

  double operator()(const std::vector<double>& x) const;

  double Up() const;

private:
  std::shared_ptr<ControlParameter> _myDataPtr;
  ParameterList& _parList;
};
  }  // namespace Minuit2

}  // namespace ROOT


#endif 
