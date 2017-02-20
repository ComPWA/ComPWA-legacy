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

#include <boost/property_tree/ptree.hpp>

#include "Core/DataPoint.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Spin.hpp"
#include "Core/FunctionTree.hpp"
//#include "Physics/HelicityFormalism/RelativisticBreitWigner.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

enum normStyle {
  none, /*!< no normaliztion between Amplitudes. */
  /*!< all amplitudes are normalized to one.
   *  The normalization factor is \f$ 1/\sqrt(\int |A|^2)\f$ */
  one
};

class
    RelativisticBreitWigner; /* Forward declaration to avoid circular includes
                                */

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

  /**
   Factory for RelativisticBreitWigner dynamical function

   @param pt Configuration tree
   @return Constructed object
   */
  static std::shared_ptr<AbstractDynamicalFunction>
  Factory(boost::property_tree::ptree &pt) {
    std::shared_ptr<AbstractDynamicalFunction> ret;
    return ret;
  }

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

  //! Resonance sub system
  unsigned int _dataIndex;

  //! Resonance spin
  ComPWA::Spin _spin;

private:
  //! Resonance shape was modified (recalculate the normalization)
  bool _modified;

  //! Integral value (temporary)
  double _integral;

  //! Temporary value of mass (used to trigger recalculation of normalization)
  double _current_mass;
};

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_ABSTRACTDYNAMICALFUNCTION_HPP_ */
