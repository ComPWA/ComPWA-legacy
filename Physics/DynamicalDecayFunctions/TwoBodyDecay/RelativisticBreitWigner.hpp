//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef PHYSICS_HELICITYAMPLITUDE_RELATIVISTICBREITWIGNER_HPP_
#define PHYSICS_HELICITYAMPLITUDE_RELATIVISTICBREITWIGNER_HPP_

#include <vector>

#include "Core/Functions.hpp"
#include "Core/Exceptions.hpp"
#include "Physics/DynamicalDecayFunctions/AbstractDynamicalFunction.hpp"

class Spin;

namespace DynamicalFunctions {

/**
 * Relativistic Breit-Wigner
 * (Breit Wigner with Blatt-Weisskopf barrier factors)
 *
 * The dynamical function implemented here is taken from PDG2014 (Eq.47-22) for
 * the one channel case.
 *
 * The three required parameters are defined in this model and can be
 * retrieved via the #getParameterList() function:
 * @resonance_width_ width of the resonance
 * @resonance_mass_ mass of the resonance
 * @meson_radius_ Scale of interaction range
 *
 * The remaining required parameter is the angular momentum @J_ of the two
 * particle state.
 *
 * Additional informations from the event are extracted from the data point
 * within the #evaluate() function. For example the invariant mass and the
 * daughter masses.
 */
class RelativisticBreitWigner: public AbstractDynamicalFunction {
public:
  RelativisticBreitWigner(const Spin& J);
  virtual ~RelativisticBreitWigner();

  void initialiseParameters(const boost::property_tree::ptree& parameter_info);

  std::complex<double> evaluate(const dataPoint& point,
      unsigned int evaluation_index) const;

private:
  std::shared_ptr<DoubleParameter> resonance_width_;
  std::shared_ptr<DoubleParameter> resonance_mass_;
  std::shared_ptr<DoubleParameter> meson_radius_;
  Spin J_;
};

}

#endif /* PHYSICS_HELICITYAMPLITUDE_RELATIVISTICBREITWIGNER_HPP_ */
