//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//----------------------------------------------------------------------------------

#include "AbstractDynamicalFunction.hpp"

namespace ComPWA {
namespace Physics {
namespace DynamicalFunctions {

AbstractDynamicalFunction::AbstractDynamicalFunction() {
  // TODO Auto-generated constructor stub

}

AbstractDynamicalFunction::~AbstractDynamicalFunction() {
  // TODO Auto-generated destructor stub
}

const ParameterList& AbstractDynamicalFunction::getParameterList() const {
  return parameter_list_;
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
