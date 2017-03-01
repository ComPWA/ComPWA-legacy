
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
  CoherentIntensity();
  virtual ~CoherentIntensity();

  virtual double GetMaximum(std::shared_ptr<ComPWA::Generator> gen) const;

  virtual double Intensity(const ComPWA::dataPoint &point) const;

  virtual double IntensityNoEff(const ComPWA::dataPoint &point) const;

  //========== FUNCTIONTREE =============
  //! Check of tree is available
  virtual bool HasTree() { return 0; }

  //! Getter function for basic amp tree
  virtual std::shared_ptr<ComPWA::FunctionTree>
  GetTree(ComPWA::ParameterList &, ComPWA::ParameterList &,
          ComPWA::ParameterList &, std::string suffix="") {
    return std::shared_ptr<ComPWA::FunctionTree>();
  }

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
  Factory(const boost::property_tree::ptree &pt) {
    auto obj = std::shared_ptr<CoherentIntensity>();
    obj->SetName(pt.get<string>("AmpIntensity.<xmlattr>.Name", "empty"));
    try {
      auto strength = ComPWA::DoubleParameterFactory(pt.get_child("Strength"));
      obj->SetStrength( std::make_shared<DoubleParameter>(strength) );
    } catch (boost::property_tree::ptree_bad_path &ex) {
      /* strength is optional */
      obj->SetStrength( std::make_shared<ComPWA::DoubleParameter>("", 1.0) );
    }

    for (const auto &v : pt.get_child("Amplitude")) {
      obj->Add(
          ComPWA::Physics::HelicityFormalism::SequentialTwoBodyDecay::Factory(
              v.second));
    }
    return obj;
  }

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
  virtual double Integral() const;


  std::vector<std::shared_ptr<
      ComPWA::Physics::HelicityFormalism::SequentialTwoBodyDecay>>
      _seqDecays;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_ */
