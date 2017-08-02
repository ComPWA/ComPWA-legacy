//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		 Peter Weidenkaff - adding UnitAmp
//-------------------------------------------------------------------------------

#ifndef AMPINTENSITY_HPP_
#define AMPINTENSITY_HPP_

#include <vector>
#include <memory>
#include <math.h>

#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "DataReader/Data.hpp"

#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Generator.hpp"

namespace ComPWA {

/*! \class AmpIntensity
 * This class provides the interface an amplitude intensity. The intensity can
 * be moduled with a (double) strength parameter.
 * Since the intensity is a physically observable quantity it has to be
 * corrected for the space depended reconstruction efficiency. The normalization
 * has to take this into account as well.
 */
class AmpIntensity {

public:
  //============ CONSTRUCTION ==================

  //! Constructor with an optional name, strength and efficiency
  AmpIntensity(std::string name = "",
               std::shared_ptr<DoubleParameter> strength =
                   std::shared_ptr<DoubleParameter>(new DoubleParameter("",
                                                                        1.0)),
               std::shared_ptr<Efficiency> eff =
                   std::shared_ptr<Efficiency>(new UnitEfficiency))
      : _name(name), _eff(eff), _strength(strength), _current_strength(-999) {
    _strength->FixParameter(true);
  }

  virtual ~AmpIntensity() { /* nothing */
  }

  //! Function to create a full copy
  virtual AmpIntensity *Clone(std::string newName = "") const = 0;

  //================ EVALUATION =================

  /*! Evaluate intensity at dataPoint in phase-space
   * @param point Data point
   * @return Intensity
   */
  virtual double Intensity(const dataPoint &point) const = 0;

  //============ SET/GET =================
  //! Get name
  virtual std::string Name() const { return _name; }

  //! Get strength parameter
  double Strength() const { return _strength->GetValue(); }

  virtual void GetParameters(ParameterList &list) = 0;
  
  //! Fill vector with parameters
  virtual void GetParametersFast(std::vector<double> &list) const {
    list.push_back(_strength->GetValue());
  }

  //! Fill ParameterList with fit fractions
  virtual void GetFitFractions(ParameterList &parList) = 0;

  /*! Set phase space samples
   * We use phase space samples to calculate the normalizations. In case of
   * intensities we phase space sample phspSample includes the event efficiency.
   * The sample toySample is used for normalization calculation for e.g.
   * Resonacnes without efficiency.
   */
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::dataPoint>> toySample) = 0;

  virtual std::shared_ptr<AmpIntensity> GetComponent(std::string name) = 0;
  
  virtual void Reset() {};
  //========== FUNCTIONTREE =============

  //! Check of tree is available
  virtual bool HasTree() const { return false; }

  /*! Get FunctionTree
   * @param sample Data sample
   * @param phspSample Sample of phase space distributed events including
   * efficiency.
   * @param toySample Sample of phase space distributed events without
   * efficiency.
   */
  virtual std::shared_ptr<FunctionTree> GetTree(const ParameterList &sample,
                                                const ParameterList &phspSample,
                                                const ParameterList &toySample,
                                                std::string suffix = "") {
    return std::shared_ptr<FunctionTree>();
  }

  //======== ITERATORS/OPERATORS =============

public:

protected:
  //! Name
  std::string _name;

  //! Phase space depended efficiency
  std::shared_ptr<Efficiency> _eff;

  std::shared_ptr<ComPWA::DoubleParameter> _strength;
  //! temporary strength
  double _current_strength;
};
//-----------------------------------------------------------------------------

//! Split string into pieces which are separated by blanks
//Todo: wohin?
static std::vector<std::string> splitString(std::string str) {
std::vector<std::string> result;
std::istringstream iStr(str);
std::vector<std::string> stringFrag{std::istream_iterator<std::string>{iStr},
                                    std::istream_iterator<std::string>{}};
for (auto i : stringFrag) {
  result.push_back(i);
}
return result;
}

} /* namespace ComPWA */
#endif
