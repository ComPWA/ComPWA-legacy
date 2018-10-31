// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/AmpIntensity.hpp"
#include "Physics/CoherentIntensity.hpp"
#include "Tools/Integration.hpp"

#ifndef INCOHERENT_INTENSITY_HPP
#define INCOHERENT_INTENSITY_HPP

namespace ComPWA {
namespace Physics {

class IncoherentIntensity : public ComPWA::AmpIntensity,
      public std::enable_shared_from_this<IncoherentIntensity> {

public:
  IncoherentIntensity() : ComPWA::AmpIntensity(), PhspVolume(1.0) {}

  IncoherentIntensity(std::shared_ptr<PartList> partL,
                      std::shared_ptr<Kinematics> kin,
                      const boost::property_tree::ptree &pt);

  /// Function to create a full copy of the amplitude
  ComPWA::AmpIntensity *clone(std::string newName = "") const {
    auto tmp = (new IncoherentIntensity(*this));
    tmp->Name = newName;
    return tmp;
  }

  void load(std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
          const boost::property_tree::ptree &pt);

  virtual boost::property_tree::ptree save() const;

  /// Calculate intensity of amplitude at point in phase-space
  virtual double intensity(const ComPWA::DataPoint &point) const;

  void addIntensity(std::shared_ptr<ComPWA::AmpIntensity> intens) {
    Intensities.push_back(intens);
  }

  /// Summand of incoherent sum.
  std::shared_ptr<ComPWA::AmpIntensity> intensity(int pos) {
    return Intensities.at(pos);
  }

  /// List of summands of incoherent sum.
  std::vector<std::shared_ptr<ComPWA::AmpIntensity>> &intensities() {
    return Intensities;
  }

  virtual void reset() {
    Intensities.clear();
//    NormalizationValues.clear();
//    Parameters.clear();
    return;
  }

  /// Add parameters to \p list only if not already in.
  virtual void parameters(ComPWA::ParameterList &list);

  /// Fill vector with parameters.
  virtual void parametersFast(std::vector<double> &list) const {
    AmpIntensity::parametersFast(list);
    for (auto i : Intensities) {
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
  virtual void
  setPhspSample(std::shared_ptr<std::vector<ComPWA::DataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::DataPoint>> toySample) {
    PhspSample = phspSample;
    for (auto i : Intensities) {
      i->setPhspSample(phspSample, toySample);
    }
  };

  virtual void setPhspVolume(double vol) { PhspVolume = vol; };

  virtual std::shared_ptr<AmpIntensity> component(std::string name);
  /// get the component of amplitude of decay resName -> daug1 Name + daug2Name 
  /// in the decay, orbitan angular momentum = L and spin = S
  /// if L/S < 0, then all possible L/S are included
  /// if daug1Name and/or daug2Name == "", then all possbile decays/decay with
  /// one daughter is daug2/daug1 are included.
  virtual std::shared_ptr<AmpIntensity> component(std::string name,
      std::string resName, std::string daug1Name = "", std::string daug2Name = "",
      int L = -1, int S = -1);


  /// Check of tree is available
  virtual bool hasTree() const { return true; }

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
  std::vector<std::vector<double>> Parameters;

  std::vector<std::shared_ptr<ComPWA::AmpIntensity>> Intensities;
};

} // namespace Physics
} // namespace ComPWA

#endif
