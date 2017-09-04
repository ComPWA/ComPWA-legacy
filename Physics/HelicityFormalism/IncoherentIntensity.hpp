// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/AmpIntensity.hpp"
#include "Physics/HelicityFormalism/CoherentIntensity.hpp"
#include "Tools/Integration.hpp"

#ifndef INCOHERENT_INTENSITY_HPP
#define INCOHERENT_INTENSITY_HPP

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class IncoherentIntensity : public ComPWA::AmpIntensity {

public:
  //============ CONSTRUCTION ==================

  IncoherentIntensity() : ComPWA::AmpIntensity(), phspVolume_(1.0) {}

  //! Function to create a full copy of the amplitude
  ComPWA::AmpIntensity *Clone(std::string newName = "") const {
    auto tmp = (new IncoherentIntensity(*this));
    tmp->_name = newName;
    return tmp;
  }

  static std::shared_ptr<IncoherentIntensity>
  Factory(std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
          const boost::property_tree::ptree &pt);

  static boost::property_tree::ptree
  Save(std::shared_ptr<IncoherentIntensity> intens);

  //================ EVALUATION =================

  /// Calculate intensity of amplitude at point in phase-space
  virtual double Intensity(const ComPWA::dataPoint &point) const;

  //================== SET/GET =================

  void AddIntensity(std::shared_ptr<ComPWA::AmpIntensity> intens) {
    _intens.push_back(intens);
  }

  std::shared_ptr<ComPWA::AmpIntensity> GetIntensity(int pos) {
    return _intens.at(pos);
  }

  std::vector<std::shared_ptr<ComPWA::AmpIntensity>> &GetIntensities() {
    return _intens;
  }

  virtual void Reset() {
    _intens.clear();
    return;
  }

  /// Add parameters to \p list.
  /// Add parameters to list only if not already in
  virtual void GetParameters(ComPWA::ParameterList &list);

  /// Fill vector with parameters
  virtual void GetParametersFast(std::vector<double> &list) const {
    AmpIntensity::GetParametersFast(list);
    for (auto i : _intens) {
      i->GetParametersFast(list);
    }
  }
  
  /// Update parameters in AmpIntensity to the values given in \p list
  virtual void UpdateParameters(const ParameterList &list);

  /// Calculate & fill fit fractions of this amplitude to ParameterList
  virtual void GetFitFractions(ComPWA::ParameterList &parList){};

  /// Set phase space sample.
  /// We use a phase space sample to calculate the normalization and determine
  /// the maximum of the amplitude. In case that the efficiency is already
  /// applied to the sample set fEff to false.
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::dataPoint>> toySample) {
    _phspSample = phspSample;
    for (auto i : _intens) {
      i->SetPhspSample(phspSample, toySample);
    }
  };

  virtual void SetPhspVolume(double vol) { phspVolume_ = vol; };

  virtual std::shared_ptr<AmpIntensity> GetComponent(std::string name);

  //======== ITERATORS/OPERATORS =============
  typedef std::vector<std::shared_ptr<ComPWA::AmpIntensity>>::iterator
      coherentIntItr;

  coherentIntItr First() { return _intens.begin(); }

  coherentIntItr Last() { return _intens.end(); }

  //=========== FUNCTIONTREE =================

  /// Check of tree is available
  virtual bool HasTree() const { return true; }

  /// Get FunctionTree
  virtual std::shared_ptr<ComPWA::FunctionTree>
  GetTree(std::shared_ptr<Kinematics> kin, const ComPWA::ParameterList &sample,
          const ComPWA::ParameterList &phspSample,
          const ComPWA::ParameterList &toySample, unsigned int nEvtVar,
          std::string suffix = "");

protected:
  /// Phase space sample to calculate the normalization and maximum value.
  std::shared_ptr<std::vector<ComPWA::dataPoint>> _phspSample;

  double phspVolume_;

  /// Caching of normalization values
  std::vector<double> _normValues;

  /// Temporary storage of the para
  std::vector<std::vector<double>> _parameters;

  std::vector<std::shared_ptr<ComPWA::AmpIntensity>> _intens;
};

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA

#endif
