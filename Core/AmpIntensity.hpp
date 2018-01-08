// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// AmpIntensity base class.
///

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

///
/// \class AmpIntensity
/// This class provides the interface an amplitude intensity. The intensity can
/// be moduled with a (double) strength parameter.
/// Since the intensity is a physically observable quantity it has to be
/// corrected for the space depended reconstruction efficiency. The
/// normalization
/// has to take this into account as well.
///
/// \par Efficiency correction
/// ToDo: explain unbinned versus binned efficiency correction
///
class AmpIntensity {

public:
  /// Constructor with an optional name, strength and efficiency
  AmpIntensity(std::string name = "",
               std::shared_ptr<FitParameter> strength =
                   std::shared_ptr<FitParameter>(new FitParameter("", 1.0)),
               std::shared_ptr<ComPWA::Efficiency> eff =
                   std::shared_ptr<ComPWA::Efficiency>(new UnitEfficiency))
      : Name(name), Eff(eff), Strength(strength), CurrentStrength(-999) {
    Strength->fixParameter(true);
  }

  virtual ~AmpIntensity(){};

  /// Clone pattern. Create a clone object with \p newName of type Base class.
  virtual AmpIntensity *clone(std::string newName = "") const = 0;

  virtual boost::property_tree::ptree save() const = 0;

  /// Evaluate intensity of model at \p point in phase-space
  virtual double intensity(const DataPoint &point) const = 0;

  virtual std::string name() const { return Name; }

  virtual void setName(std::string n) { Name = n; }

  /// Get strength parameter
  double strength() const { return Strength->value(); }

  virtual void parameters(ParameterList &list) = 0;

  /// Fill vector with parameters
  virtual void parametersFast(std::vector<double> &list) const {
    list.push_back(Strength->value());
  }

  /// Update parameters in AmpIntensity to the values given in \p list
  virtual void updateParameters(const ParameterList &list) = 0;

  /// Set phase space samples
  /// We use phase space samples to calculate the normalizations. In case of
  /// intensities we phase space sample phspSample includes the event
  /// efficiency.
  /// The sample toySample is used for normalization calculation for e.g.
  /// Resonacnes without efficiency.
  virtual void
  setPhspSample(std::shared_ptr<std::vector<ComPWA::DataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::DataPoint>> toySample) = 0;

  virtual std::shared_ptr<AmpIntensity> component(std::string name) = 0;

  virtual void reset(){};
  //========== FUNCTIONTREE =============

  /// Check if a FunctionTree is implemented for a certain (derived) class.
  virtual bool hasTree() const { return false; }

  virtual std::shared_ptr<FunctionTree>
  tree(std::shared_ptr<Kinematics> kin, const ParameterList &sample,
       const ParameterList &phspSample, const ParameterList &toySample,
       unsigned int nEvtVar, std::string suffix = "") = 0;

protected:
  /// Name
  std::string Name;

  /// Phase space depended efficiency
  std::shared_ptr<Efficiency> Eff;

  std::shared_ptr<ComPWA::FitParameter> Strength;
  //    std::shared_ptr<ComPWA::FunctionTree> Strength;

  /// temporary strength
  double CurrentStrength;
};

/// Split string into pieces which are separated by blanks
inline std::vector<std::string> splitString(std::string str) {
  std::vector<std::string> result;
  std::istringstream iStr(str);
  std::vector<std::string> stringFrag{std::istream_iterator<std::string>{iStr},
                                      std::istream_iterator<std::string>{}};
  for (auto i : stringFrag) {
    result.push_back(i);
  }
  return result;
}

} // ns::ComPWA
#endif
