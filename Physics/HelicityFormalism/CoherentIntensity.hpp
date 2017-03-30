
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
#include "Physics/HelicityFormalism/SequentialTwoBodyDecay.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class CoherentIntensity : public ComPWA::AmpIntensity {

public:
  CoherentIntensity(){};
  virtual ~CoherentIntensity(){};

  virtual double GetMaximum(std::shared_ptr<ComPWA::Generator> gen) const {
    return 1.0;
  }

  virtual double Intensity(const ComPWA::dataPoint &point) const;

  virtual double IntensityNoEff(const ComPWA::dataPoint &point) const;

  //========== FUNCTIONTREE =============
  //! Check of tree is available
  virtual bool HasTree() { return 0; }

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

  void Add(std::shared_ptr<
           ComPWA::Physics::HelicityFormalism::SequentialTwoBodyDecay>
               d) {
    _seqDecays.push_back(d);
  }

  //! Get number of partial decays
  size_t size() const { return _seqDecays.size(); }

  typedef std::vector<std::shared_ptr<
      ComPWA::Physics::HelicityFormalism::SequentialTwoBodyDecay>>::iterator
      seqDecayItr;

  seqDecayItr First() { return _seqDecays.begin(); }

  seqDecayItr Last() { return _seqDecays.end(); }

  //=========== PARAMETERS =================
  /**! Add amplitude parameters to list
   * Add parameters only to list if not already in
   * @param list Parameter list to be filled
   */
  virtual void GetParameters(ComPWA::ParameterList &list) const {};

  //! Calculate & fill fit fractions of this amplitude to ParameterList
  virtual void GetFitFractions(ComPWA::ParameterList &parList){};

  //============= PRINTING =====================
  //! Print amplitude to logging system
  virtual void to_str() const { LOG(info) << "CoherentIntensity"; }

protected:
  virtual double Integral() const { return 1.0; }

  virtual std::shared_ptr<FunctionTree>
  setupBasicTree(ParameterList &sample, ParameterList &phspSample,
                 std::string suffix="") const;

  std::vector<std::shared_ptr<
      ComPWA::Physics::HelicityFormalism::SequentialTwoBodyDecay>>
      _seqDecays;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_ */
