// Copyright (c) 2018 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/AmpIntensity.hpp"
#include "Core/Logging.hpp"
#include "Physics/EvtGen/DalitzKinematics.hpp"
#include "Physics/EvtGen/EvtDalitzPlot.hh"
#include "Physics/EvtGen/EvtDalitzReso.hh"
#include "Physics/SubSystem.hpp"
#include "Tools/Integration.hpp"

#ifndef COMPWAEVTGENIF_HPP
#define COMPWAEVTGENIF_HPP

namespace ComPWA {
namespace Physics {
namespace EvtGenIF {

class EvtGenIF : public ComPWA::AmpIntensity {

public:
  EvtGenIF() : ComPWA::AmpIntensity(), PhspVolume(1.0), DalitzPlot() {}

  EvtGenIF(std::shared_ptr<PartList> partL, double mA, double mB, double mC,
           double bigM, double ldel = 0., double rdel = 0.)
      : ComPWA::AmpIntensity(), PhspVolume(1.0),
        DalitzPlot(mA, mB, mC, bigM, ldel, rdel) {}

  /// Function to create a full copy of the amplitude
  ComPWA::AmpIntensity *clone(std::string newName = "") const {
    auto tmp = (new EvtGenIF(*this));
    tmp->Name = newName;
    return tmp;
  }

  virtual boost::property_tree::ptree save() const;

  virtual std::shared_ptr<AmpIntensity> component(std::string name);
  virtual std::shared_ptr<AmpIntensity> component(std::string name,
      std::string resName, std::string daug1Name, std::string daug2Name,
      int L, int S);

  /// Add EvtGen Dalitz Resonance
  virtual void addResonance(std::string name, double m0, double g0, double spin,
                            ComPWA::Physics::SubSystem subsys);

  /// Add EvtGen Dalitz Resonance
  virtual void addHeliResonance(boost::property_tree::ptree pt,
                                std::shared_ptr<PartList> partL);

  /// Add EvtGen Dalitz Resonances from XML model
  virtual void addResonances(boost::property_tree::ptree pt,
                             std::shared_ptr<DalitzKinematics> kin,
                             std::shared_ptr<PartList> partL);

  /// Calculate intensity of amplitude at point in phase-space
  virtual double intensity(const ComPWA::DataPoint &point) const;

  virtual void reset() {
    //    NormalizationValues.clear();
    //    Parameters.clear();
    return;
  }

  /// Add parameters to \p list only if not already in.
  virtual void parameters(ComPWA::ParameterList &list);

  /// Fill vector with parameters.
  virtual void parametersFast(std::vector<double> &list) const {
    AmpIntensity::parametersFast(list);
  }

  /// Update parameters in AmpIntensity to the values given in \p list
  virtual void updateParameters(const ParameterList &list);

  /// Set phase space sample.
  /// We use a phase space sample to calculate the normalization and determine
  /// the maximum of the amplitude. In case that the efficiency is already
  /// applied to the sample set fEff to false.
  virtual void
  setPhspSample(std::shared_ptr<std::vector<ComPWA::DataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::DataPoint>> toySample) {
    PhspSample = phspSample;
  };

  virtual void setPhspVolume(double vol) { PhspVolume = vol; };

  /// Check of tree is available
  virtual bool hasTree() const { return false; }

  /// Get FunctionTree
  virtual std::shared_ptr<ComPWA::FunctionTree>
  tree(std::shared_ptr<Kinematics> kin, const ComPWA::ParameterList &sample,
       const ComPWA::ParameterList &phspSample,
       const ComPWA::ParameterList &toySample, unsigned int nEvtVar,
       std::string suffix = "");

protected:
  /// Phase space sample to calculate the normalization and maximum value.
  std::shared_ptr<std::vector<ComPWA::DataPoint>> PhspSample;

  double PhspVolume;

  /// Caching of normalization values for each intensity
  std::vector<double> NormalizationValues;

  /// Temporary storage of the para
  // std::vector<std::vector<double>> Parameters;
  std::map<std::string, std::shared_ptr<ComPWA::FitParameter>> evtPars;

  EvtDalitzPlot DalitzPlot;
  std::vector<EvtDalitzReso> Resos;

  EvtDalitzPoint transformToEvt(const ComPWA::DataPoint &point) const {
    return EvtDalitzPoint(point.value(0), point.value(1), point.value(2),
                          point.value(3), point.value(4), point.value(5));
  }
};

} // namespace EvtGenIF
} // namespace Physics
} // namespace ComPWA

#endif
