// Copyright (c) 2018 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_EVTGEN_EVTGENIF_HPP
#define COMPWA_PHYSICS_EVTGEN_EVTGENIF_HPP

#include "Core/Event.hpp"
#include "Core/FunctionTree/FitParameter.hpp"
#include "Core/Logging.hpp"
#include "Data/DataSet.hpp"
#include "Physics/EvtGen/DalitzKinematics.hpp"
#include "Physics/SubSystem.hpp"
#include "ThirdParty/EvtGen/EvtDalitzPlot.hh"
#include "ThirdParty/EvtGen/EvtDalitzReso.hh"
#include "Tools/Integration.hpp"

namespace ComPWA {
namespace Physics {
namespace EvtGen {

class EvtGenIF : public ComPWA::Intensity {

public:
  EvtGenIF() : ComPWA::Intensity(), PhspVolume(1.0), DalitzPlot() {}

  EvtGenIF(double mA, double mB, double mC, double bigM, double ldel = 0.,
           double rdel = 0.)
      : ComPWA::Intensity(), PhspVolume(1.0),
        DalitzPlot(mA, mB, mC, bigM, ldel, rdel) {}

  /// Add EvtGen Dalitz Resonance
  void addResonance(const std::string &name, double m0, double g0, double spin,
                    const ComPWA::Physics::SubSystem &subsys);

  /// Add EvtGen Dalitz Resonance
  void addHeliResonance(const boost::property_tree::ptree &pt,
                        const ComPWA::ParticleList &partL);

  /// Add EvtGen Dalitz Resonances from XML model
  void addResonances(const boost::property_tree::ptree &pt,
                     std::shared_ptr<DalitzKinematics> kin,
                     const ComPWA::ParticleList &partL);

  std::vector<double> evaluate(const ComPWA::DataMap &data) noexcept;

  /// It is important to input the vector in the same length and order as
  /// defined in the getParameters() method. So in other words, call
  /// getParameters() first, then modify the contents and finally input them in
  /// this method.
  void updateParametersFrom(const std::vector<double> &Parameters) final;
  std::vector<ComPWA::Parameter> getParameters() const final;

  /// Set phase space sample.
  /// We use a phase space sample to calculate the normalization and determine
  /// the maximum of the amplitude. In case that the efficiency is already
  /// applied to the sample set fEff to false.
  virtual void setPhspSample(std::shared_ptr<ComPWA::Data::DataSet> phspSample,
                             std::shared_ptr<ComPWA::Data::DataSet> toySample) {
    PhspSample = phspSample;
  };

  virtual void setPhspVolume(double vol) { PhspVolume = vol; };

private:
  /// Phase space sample to calculate the normalization and maximum value.
  std::shared_ptr<ComPWA::Data::DataSet> PhspSample;

  double PhspVolume;

  /// Caching of normalization values for each intensity
  std::vector<double> NormalizationValues;

  /// Temporary storage of the para
  // std::vector<std::vector<double>> Parameters;
  std::map<std::string, std::shared_ptr<ComPWA::FunctionTree::FitParameter>>
      evtPars;

  EvtDalitzPlot DalitzPlot;
  std::vector<EvtDalitzReso> Resos;
};

} // namespace EvtGen
} // namespace Physics
} // namespace ComPWA

#endif
