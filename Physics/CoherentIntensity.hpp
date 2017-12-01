// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_
#define PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/AmpIntensity.hpp"
#include "DataReader/Data.hpp"
#include "Physics/SequentialPartialAmplitude.hpp"

namespace ComPWA {
namespace Physics {

class CoherentIntensity : public ComPWA::AmpIntensity {

public:
  //============ CONSTRUCTION ==================

  CoherentIntensity(
      std::string name = "",
      std::shared_ptr<DoubleParameter> strength =
          std::shared_ptr<DoubleParameter>(new DoubleParameter("", 1.0)),
      std::shared_ptr<Efficiency> eff =
          std::shared_ptr<Efficiency>(new UnitEfficiency))
      : AmpIntensity(name, strength, eff), PhspVolume(1.0){};

  virtual ~CoherentIntensity(){};

  /// Clone pattern
  ComPWA::AmpIntensity *clone(std::string newName = "") const {
    auto tmp = (new CoherentIntensity(*this));
    tmp->Name = newName;
    return tmp;
  }

  static std::shared_ptr<CoherentIntensity>
  Factory(std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
          const boost::property_tree::ptree &pt);

  static boost::property_tree::ptree
  Save(std::shared_ptr<CoherentIntensity> intens);

  //================ EVALUATION =================

  /// Calculate intensity of amplitude at point in phase-space
  virtual double intensity(const ComPWA::DataPoint &point) const;

  //============ SET/GET =================

  void addAmplitude(std::shared_ptr<ComPWA::Physics::Amplitude> decay) {
    Amplitudes.push_back(decay);
  }

  std::shared_ptr<ComPWA::Physics::Amplitude> amplitude(int pos) {
    return Amplitudes.at(pos);
  }

  std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>> &amplitudes() {
    return Amplitudes;
  }

  virtual void reset() {
    Amplitudes.clear();
    return;
  }

  /// Add parameters to \p list.
  /// Add parameters to list only if not already in
  virtual void parameters(ComPWA::ParameterList &list);

  /// Fill vector with parameters
  virtual void parametersFast(std::vector<double> &list) const {
    AmpIntensity::parametersFast(list);
    for (auto i : Amplitudes) {
      i->parametersFast(list);
    }
  }
  
  /// Update parameters in AmpIntensity to the values given in \p list
  virtual void updateParameters(const ParameterList &list);

  /// Calculate & fill fit fractions of this amplitude to ParameterList
  virtual void fitFractions(ComPWA::ParameterList &parList){};

  /// Set phase space sample.
  /// We use a phase space sample to calculate the normalization and determine
  /// the maximum of the amplitude. In case that the efficiency is already
  /// applied to the sample set fEff to false.
  void
  setPhspSample(std::shared_ptr<std::vector<ComPWA::DataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::DataPoint>> toySample) {
    PhspSample = phspSample;

    for (auto i : Amplitudes)
      i->setPhspSample(toySample);
  };
  
  virtual void setPhspVolume(double vol) { PhspVolume = vol; };

  virtual std::shared_ptr<AmpIntensity> component(std::string name);

  //======== ITERATORS/OPERATORS =============

  typedef std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>>::iterator
      seqDecayItr;

  seqDecayItr first() { return Amplitudes.begin(); }

  seqDecayItr last() { return Amplitudes.end(); }

  //========== FUNCTIONTREE =============

  /// Check of tree is available
  virtual bool hasTree() const { return true; }

  /// Getter function for basic amp tree
  virtual std::shared_ptr<ComPWA::FunctionTree>
  tree(std::shared_ptr<Kinematics> kin, const ComPWA::ParameterList &sample,
          const ComPWA::ParameterList &phspSample,
          const ComPWA::ParameterList &toySample, unsigned int nEvtVar,
          std::string suffix = "");

protected:
  /// Phase space sample to calculate the normalization and maximum value.
  std::shared_ptr<std::vector<ComPWA::DataPoint>> PhspSample;
  
  double PhspVolume;

  virtual std::shared_ptr<FunctionTree>
  setupBasicTree(std::shared_ptr<Kinematics> kin, const ParameterList &sample,
                 const ParameterList &phspSample,
                 std::string suffix = "") const;

  std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>> Amplitudes;

};

} // namespace Physics
} // namespace ComPWA

#endif
