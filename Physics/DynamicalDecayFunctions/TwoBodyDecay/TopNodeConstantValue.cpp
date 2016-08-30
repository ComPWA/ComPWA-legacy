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

#include "Physics/DynamicalDecayFunctions/TwoBodyDecay/TopNodeConstantValue.hpp"

namespace ComPWA {
namespace Physics {
namespace DynamicalFunctions {

TopNodeConstantValue::TopNodeConstantValue() {
  // TODO Auto-generated constructor stub

}

TopNodeConstantValue::~TopNodeConstantValue() {
  // TODO Auto-generated destructor stub
}

void TopNodeConstantValue::initialiseParameters(
    const boost::property_tree::ptree& parameter_info,
    const ExternalParameters& external_parameters) {
}

std::complex<double> TopNodeConstantValue::evaluate(const dataPoint& point,
    unsigned int evaluation_index) const {
  return std::complex<double>(1.0, 0.0);
}

std::complex<double> TopNodeConstantValue::evaluate(unsigned int storage_index,
    unsigned int data_index, unsigned int evaluation_index) const {
  return std::complex<double>(1.0, 0.0);
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
