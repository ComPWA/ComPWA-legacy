// Copyright (c) 2018 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Event.hpp"
#include "Core/FunctionTree/Intensity.hpp"
#include "Core/Logging.hpp"
#include "Physics/EvtGen/DalitzKinematics.hpp"
#include "Physics/SubSystem.hpp"
#include "ThirdParty/EvtGen/EvtDalitzPlot.hh"
#include "ThirdParty/EvtGen/EvtDalitzReso.hh"
#include "Tools/Integration.hpp"

#ifndef COMPWA_PHYSICS_EVTGEN_EVTGENIF_HPP
#define COMPWA_PHYSICS_EVTGEN_EVTGENIF_HPP

namespace ComPWA {
namespace Physics {
namespace EvtGen {

class EvtGenIF : public ComPWA::Intensity {

public:
  EvtGenIF() : ComPWA::Intensity(), PhspVolume(1.0), DalitzPlot() {}

  EvtGenIF(std::shared_ptr<PartList> partL, double mA, double mB, double mC,
           double bigM, double ldel = 0., double rdel = 0.)
      : ComPWA::Intensity(), PhspVolume(1.0),
        DalitzPlot(mA, mB, mC, bigM, ldel, rdel) {}

  virtual boost::property_tree::ptree save() const;

  virtual std::shared_ptr<Intensity> component(const std::string &name);

  /// Add EvtGen Dalitz Resonance
  virtual void addResonance(const std::string &name, double m0, double g0,
                            double spin,
                            const ComPWA::Physics::SubSystem &subsys);

  /// Add EvtGen Dalitz Resonance
  virtual void addHeliResonance(const boost::property_tree::ptree &pt,
                                std::shared_ptr<PartList> partL);

  /// Add EvtGen Dalitz Resonances from XML model
  virtual void addResonances(const boost::property_tree::ptree &pt,
                             std::shared_ptr<DalitzKinematics> kin,
                             std::shared_ptr<PartList> partL);

  /// Calculate intensity of amplitude at point in phase-space
  virtual double evaluate(const ComPWA::DataPoint &point) const;

  virtual void reset() {
    //    NormalizationValues.clear();
    //    Parameters.clear();
    return;
  }

  /// Add parameters to \p list only if not already in.
  virtual void addUniqueParametersTo(ComPWA::FunctionTree::ParameterList &list);

  /// Update parameters in AmpIntensity to the values given in \p list
  virtual void
  updateParametersFrom(const ComPWA::FunctionTree::ParameterList &list);

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

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createFunctionTree(const ComPWA::FunctionTree::ParameterList &DataSample,
                     const std::string &suffix) const;

protected:
  /// Phase space sample to calculate the normalization and maximum value.
  std::shared_ptr<std::vector<ComPWA::DataPoint>> PhspSample;

  double PhspVolume;

  /// Caching of normalization values for each intensity
  std::vector<double> NormalizationValues;

  /// Temporary storage of the para
  // std::vector<std::vector<double>> Parameters;
  std::map<std::string, std::shared_ptr<ComPWA::FunctionTree::FitParameter>>
      evtPars;

  EvtDalitzPlot DalitzPlot;
  std::vector<EvtDalitzReso> Resos;

  EvtDalitzPoint transformToEvt(const ComPWA::DataPoint &point) const {
    return EvtDalitzPoint(
        point.KinematicVariableList[0], point.KinematicVariableList[1],
        point.KinematicVariableList[2], point.KinematicVariableList[3],
        point.KinematicVariableList[4], point.KinematicVariableList[5]);
  }
};

} // namespace EvtGen
} // namespace Physics
} // namespace ComPWA

#endif
