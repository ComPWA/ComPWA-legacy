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

#ifndef PHYSICS_HELICITYAMPLITUDE_ABSTRACTDYNAMICALFUNCTION_HPP_
#define PHYSICS_HELICITYAMPLITUDE_ABSTRACTDYNAMICALFUNCTION_HPP_

#include <complex>

#include "Core/DataPoint.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Spin.hpp"

namespace ComPWA {
namespace Physics {
namespace DynamicalFunctions {

class AbstractDynamicalFunction {
public:
  AbstractDynamicalFunction();

  virtual ~AbstractDynamicalFunction();

  virtual std::complex<double> Evaluate(const dataPoint &) const = 0;

  /**! Get current normalization.  */
  virtual double GetNormalization() = 0;
  
  /**! Setup function tree */
  virtual std::shared_ptr<FunctionTree> SetupTree(ParameterList &sample,
                                                  ParameterList &toySample,
                                                  std::string suffix) = 0;

protected:
  //! Name of resonance
  std::string _name;

  unsigned int dataIndex;

  //! Type of resonance normalization
  normStyle _normStyle;

  //! Precision of MC integration
  int _mcPrecision;

  //! Integral
  virtual double integral() { return 1.0; };

  //! Masses of daughter particles
  double _mass1, _mass2;

  //! Resonance mass
  std::shared_ptr<DoubleParameter> _mass;

  //! Barrier radi for resonance and mother particle
  std::shared_ptr<DoubleParameter> _mesonRadius;

  //! Resonance sub system
  unsigned int _dataIndex;

  //! Resonance spin
  ComPWA::Spin _spin;

private:
  //! Resonance shape was modified (recalculate the normalization)
  bool _modified;
  //! Integral value (temporary)
  double _integral;
  double _current_mass;
  double _current_mesonRadius;
};

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_ABSTRACTDYNAMICALFUNCTION_HPP_ */

/
