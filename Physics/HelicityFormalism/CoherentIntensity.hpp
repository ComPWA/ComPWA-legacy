
//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

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
  CoherentIntensity() : _maxIntens(0.), _integral(0.){};
  virtual ~CoherentIntensity(){};

  virtual double GetMaximum(std::shared_ptr<ComPWA::Generator> gen) const {
    if (_maxIntens <= 0)
      Integral(); //_maxIntens is updated in Integral()
    return _maxIntens;
  }
  //! Calculate normalization
  virtual double GetNormalization() const;

  virtual double Intensity(const ComPWA::dataPoint &point) const;

  virtual double IntensityNoEff(const ComPWA::dataPoint &point) const;

  //========== FUNCTIONTREE =============
  //! Check of tree is available
  virtual bool HasTree() const { return true; }

  //! Getter function for basic amp tree
  virtual std::shared_ptr<ComPWA::FunctionTree>
  GetTree(ComPWA::ParameterList &sample, ComPWA::ParameterList &phspSample,
          ComPWA::ParameterList &toySample, std::string suffix = "");

  /**
   Clone function

   @param newName New name
   @return Pointer to new object
   */
  ComPWA::AmpIntensity *Clone(std::string newName = "") const {
    auto tmp = (new CoherentIntensity(*this));
    tmp->SetName(newName);
    return tmp;
  }

  static std::shared_ptr<CoherentIntensity>
  Factory(const boost::property_tree::ptree &pt);

  static boost::property_tree::ptree
  Save(std::shared_ptr<CoherentIntensity> intens);

  void Add(std::shared_ptr<ComPWA::Physics::Amplitude> d) {
    _seqDecays.push_back(d);
  }

  std::shared_ptr<ComPWA::Physics::Amplitude> GetDecay(int pos) {
    return _seqDecays.at(pos);
  }

  std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>> &GetDecays() {
    return _seqDecays;
  }

  //! Get number of partial decays
  size_t size() const { return _seqDecays.size(); }

  typedef std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>>::iterator
      seqDecayItr;

  seqDecayItr First() { return _seqDecays.begin(); }

  seqDecayItr Last() { return _seqDecays.end(); }

  //=========== PARAMETERS =================
  /**! Add amplitude parameters to list
   * Add parameters only to list if not already in
   * @param list Parameter list to be filled
   */
  virtual void GetParameters(ComPWA::ParameterList &list);

  //! Calculate & fill fit fractions of this amplitude to ParameterList
  virtual void GetFitFractions(ComPWA::ParameterList &parList){};

  //============= PRINTING =====================
  //! Print amplitude to logging system
  virtual void to_str() const { LOG(info) << "CoherentIntensity"; }

  /*! Set phase space sample
   * We use a phase space sample to calculate the normalization and determine
   * the maximum of the amplitude. In case that the efficiency is already
   * applied
   * to the sample set fEff to false.
   */
  void SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample,
                     std::shared_ptr<std::vector<ComPWA::dataPoint>> toySample) {
    _phspSample = phspSample;
    
    for( auto i: _seqDecays )
      i->SetPhspSample(toySample);
  };

protected:
  virtual double Integral() const;

  virtual std::shared_ptr<FunctionTree>
  setupBasicTree(ParameterList &sample, ParameterList &phspSample,
                 std::string suffix = "") const;

  std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>> _seqDecays;

  //! Phase space sample to calculate the normalization and maximum value.
  std::shared_ptr<std::vector<ComPWA::dataPoint>> _phspSample;

  //! Maximum value
  double _maxIntens;
  double _integral;

private:
  /*! List with all parameters of the intensity.
   * We use it to check if parameters were modified and if we have to
   * recalculated the normalization.
   */
  ParameterList _currentParList;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_ */
