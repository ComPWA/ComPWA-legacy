//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#include "DynamicalFunctionFactory.hpp"

namespace DynamicalFunctions {

DynamicalFunctionFactory::DynamicalFunctionFactory() {
  // TODO Auto-generated constructor stub

}

DynamicalFunctionFactory::~DynamicalFunctionFactory() {
  // TODO Auto-generated destructor stub
}

std::shared_ptr<AbstractDynamicalFunction> DynamicalFunctionFactory::generateDynamicalFunction(
    const HelicityFormalism::TwoBodyDecayInformation& state_info) const {

}

} /* namespace DynamicalFunctions */
