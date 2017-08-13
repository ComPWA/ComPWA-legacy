// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_
#define PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/AmpIntensity.hpp"
#include "DataReader/Data.hpp"
#include "Physics/HelicityFormalism/SequentialTwoBodyDecay.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class CoherentIntensity : public ComPWA::AmpIntensity {

public:
  //============ CONSTRUCTION ==================

  CoherentIntensity(
      std::string name = "",
      std::shared_ptr<DoubleParameter> strength =
          std::shared_ptr<DoubleParameter>(new DoubleParameter("", 1.0)),
      std::shared_ptr<Efficiency> eff =
          std::shared_ptr<Efficiency>(new UnitEfficiency))
      : AmpIntensity(name, strength, eff){};

  virtual ~CoherentIntensity(){};

  //! Clone function
  ComPWA::AmpIntensity *Clone(std::string newName = "") const {
    auto tmp = (new CoherentIntensity(*this));
    tmp->_name = newName;
    return tmp;
  }

  static std::shared_ptr<CoherentIntensity>
  Factory(std::shared_ptr<PartList> partL, std::shared_ptr<Kinematics> kin,
          const boost::property_tree::ptree &pt);

  static boost::property_tree::ptree
  Save(std::shared_ptr<CoherentIntensity> intens);

  //================ EVALUATION =================

  /** Calculate intensity of amplitude at point in phase-space
   *
   * @param point Data point
   * @return
   */
  virtual double Intensity(const ComPWA::dataPoint &point) const;

  //============ SET/GET =================

  void AddAmplitude(std::shared_ptr<ComPWA::Physics::Amplitude> decay) {
    _seqDecays.push_back(decay);
  }

  std::shared_ptr<ComPWA::Physics::Amplitude> GetAmplitude(int pos) {
    return _seqDecays.at(pos);
  }

  std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>> &GetAmplitudes() {
    return _seqDecays;
  }

  virtual void Reset() {
    _seqDecays.clear();
    return;
  }

  /**! Add amplitude parameters to list
   * Add parameters only to list if not already in
   * @param list Parameter list to be filled
   */
  virtual void GetParameters(ComPWA::ParameterList &list);

  //! Fill vector with parameters
  virtual void GetParametersFast(std::vector<double> &list) const {
    AmpIntensity::GetParametersFast(list);
    for (auto i : _seqDecays) {
      i->GetParametersFast(list);
    }
  }

  //! Calculate & fill fit fractions of this amplitude to ParameterList
  virtual void GetFitFractions(ComPWA::ParameterList &parList){};

  /*! Set phase space sample
   * We use a phase space sample to calculate the normalization and determine
   * the maximum of the amplitude. In case that the efficiency is already
   * applied
   * to the sample set fEff to false.
   */
  void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::dataPoint>> toySample) {
    _phspSample = phspSample;

    for (auto i : _seqDecays)
      i->SetPhspSample(toySample);
  };
  
  virtual void SetPhspVolume(double vol) { phspVolume_ = vol; };

  virtual std::shared_ptr<AmpIntensity> GetComponent(std::string name);

  //======== ITERATORS/OPERATORS =============

  typedef std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>>::iterator
      seqDecayItr;

  seqDecayItr First() { return _seqDecays.begin(); }

  seqDecayItr Last() { return _seqDecays.end(); }

  //========== FUNCTIONTREE =============

  //! Check of tree is available
  virtual bool HasTree() const { return true; }

  //! Getter function for basic amp tree
  virtual std::shared_ptr<ComPWA::FunctionTree>
  GetTree(std::shared_ptr<Kinematics> kin, const ComPWA::ParameterList &sample,
          const ComPWA::ParameterList &phspSample,
          const ComPWA::ParameterList &toySample, unsigned int nEvtVar,
          std::string suffix = "");

protected:
  //! Phase space sample to calculate the normalization and maximum value.
  std::shared_ptr<std::vector<ComPWA::dataPoint>> _phspSample;
  
  double phspVolume_;

  virtual std::shared_ptr<FunctionTree>
  setupBasicTree(std::shared_ptr<Kinematics> kin, const ParameterList &sample,
                 const ParameterList &phspSample,
                 std::string suffix = "") const;

  std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>> _seqDecays;

};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_ */
