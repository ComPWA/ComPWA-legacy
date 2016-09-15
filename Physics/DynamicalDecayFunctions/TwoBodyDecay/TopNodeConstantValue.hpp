//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef PHYSICS_DYNAMICALDECAYFUNCTIONS_TWOBODYDECAY_TOPNODECONSTANTVALUE_HPP_
#define PHYSICS_DYNAMICALDECAYFUNCTIONS_TWOBODYDECAY_TOPNODECONSTANTVALUE_HPP_

#include "Physics/DynamicalDecayFunctions/AbstractDynamicalFunction.hpp"

namespace ComPWA {
namespace Physics {
namespace DynamicalFunctions {

class TopNodeConstantValue: public AbstractDynamicalFunction {
  void initialiseParameters(const boost::property_tree::ptree& parameter_info,
      const ExternalParameters& external_parameters);

public:
  TopNodeConstantValue();
  virtual ~TopNodeConstantValue();

  std::complex<double> evaluate(const dataPoint& point,
      unsigned int evaluation_index) const;

  std::complex<double> evaluate(unsigned int storage_index,
      unsigned int data_index, unsigned int evaluation_index) const;
};

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_DYNAMICALDECAYFUNCTIONS_TWOBODYDECAY_TOPNODECONSTANTVALUE_HPP_ */
