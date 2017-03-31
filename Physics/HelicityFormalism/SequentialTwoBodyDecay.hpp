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

#ifndef SequentialTwoBodyDecay_h
#define SequentialTwoBodyDecay_h

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/Parameter.hpp"
#include "Physics/Amplitude.hpp"
#include "Physics/HelicityFormalism/PartialDecay.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class SequentialTwoBodyDecay : Amplitude {

public:
  /**! Evaluate decay */
  virtual std::complex<double> Evaluate(const dataPoint &point) const {
    std::complex<double> result =
        std::polar(_magnitude->GetValue(), _phase->GetValue());
    for (auto i : _partDecays)
      result *= i->Evaluate(point);
    return result;
  };

  /**! Setup function tree */
  virtual std::shared_ptr<FunctionTree> GetTree(ParameterList &sample,
                                                ParameterList &phspSample,
                                                ParameterList &toySample,
                                                std::string suffix);

  /**
   Factory for SequentialTwoBodyDecay

   @param pt Configuration tree
   @return Constructed object
   */
  static std::shared_ptr<SequentialTwoBodyDecay>
  Factory(const boost::property_tree::ptree &pt);

  /**
   Add a partial decay to Sequential decay

   @param d Partial decay
   */
  void
  Add(std::shared_ptr<ComPWA::Physics::HelicityFormalism::PartialDecay> d) {
    _partDecays.push_back(d);
  }

  /**
   Get Magnitude parameter

   @return Magnitude parameter
   */
  std::shared_ptr<ComPWA::DoubleParameter> GetMagnitude() { return _magnitude; }

  /**
   Get Magnitude parameter

   @return Magnitude parameter
   */
  double GetMagnitudeValue() const { return _magnitude->GetValue(); }

  /**
   Set Magnitude parameter

   @param par Magnitude parameter
   */
  void SetMagnitude(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _magnitude = par;
  }

  /**
   Set Magnitude parameter

   @param par Magnitude parameter
   */
  void SetMagnitude(double par) { _magnitude->SetValue(par); }

  /**
   Get phase parameter

   @return Phase parameter
   */
  std::shared_ptr<ComPWA::DoubleParameter> GetPhase() { return _phase; }

  /**
   Get phase parameter

   @return Phase parameter
   */
  double GetPhaseValue() const { return _phase->GetValue(); }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  void SetPhase(std::shared_ptr<ComPWA::DoubleParameter> par) { _phase = par; }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  void SetPhase(double par) { _phase->SetValue(par); }

  //! Function to create a full copy of the amplitude
  virtual Amplitude *Clone(std::string newName = "") const {
    auto tmp = (new SequentialTwoBodyDecay(*this));
    tmp->SetName(newName);
    return tmp;
  };

  //! Print amplitude to logging system
  virtual void to_str() { LOG(info) << "SequentialTwoBodyDecay "; }

  //! Fill ParameterList with fit fractions
  virtual void GetFitFractions(ParameterList &parList){};

  /**
   Get number of partial decays

   @return Number of partial decays
   */
  size_t size() { return _partDecays.size(); };

  typedef std::vector<std::shared_ptr<PartialDecay>>::iterator partDecayItr;

  partDecayItr begin() { return _partDecays.begin(); }

  partDecayItr end() { return _partDecays.end(); }

protected:
  // TODO: we add this particle state info for the coherent sum stuff
  // the whole design is fucked up because of that, change that someday
  std::pair<ParticleStateInfo, std::pair<ParticleStateInfo, ParticleStateInfo>>
      decay_spin_info_;

  std::vector<std::shared_ptr<PartialDecay>> _partDecays;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* SequentialTwoBodyDecay_h */
