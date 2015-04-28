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

#ifndef ABSTRACTDYNAMICALFUNCTION_HPP_
#define ABSTRACTDYNAMICALFUNCTION_HPP_

#include <complex>

namespace HelicityFormalism {

class AbstractDynamicalFunction {
public:
  AbstractDynamicalFunction();
  virtual ~AbstractDynamicalFunction();

  virtual std::complex<double> evaluate() const =0;
};

} /* namespace HelicityFormalism */

#endif /* ABSTRACTDYNAMICALFUNCTION_HPP_ */
