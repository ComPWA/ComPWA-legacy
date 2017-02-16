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

#include "Core/AmpIntensity.hpp"
#include "Physics/HelicityAmplitude/SequentialTwoBodyDecay.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class CoherentIntensity : public AmpIntensity {

public:
  CoherentIntensity();
  virtual ~CoherentAmplitude();

  double GetMaximum(std::shared_ptr<Generator> gen);

  const double Intensity(const dataPoint &point);

  const double IntensityNoEff(const dataPoint &point);

  //========== FUNCTIONTREE =============
  //! Check of tree is available
  virtual bool hasTree() { return 0; }

  //! Getter function for basic amp tree
  virtual std::shared_ptr<FunctionTree>
  SetupTree(ParameterList &, ParameterList &, ParameterList &) {
    return std::shared_ptr<FunctionTree>();
  }

  AmpIntensity *Clone();

  static std::shared_ptr<CoherentIntensity>
  Factory(boost::property_tree::ptree &pt) {
    auto obj = std::make_shared<CoherentIntensity>();
    obj->SetName( pt.get<string>("AmpIntensity.<xmlattr>.name","empty") );
    BOOST_FOREACH (ptree::value_type const &ampTree,
                   pt.second.get_child("Amplitude")) {
      obj->Add( SequentialTwoBodyDecay::Factory(ampTree) );
    }
      return obj;
}
  
  void Add( std::shared_ptr<SequentialTwoBodyDecay> d ){
    _seqDecays.push_back(d);
  }

  /**
   Get number of partial decays
   
   @return Number of partial decays
   */
  size_t size() { return _seqDecayAmps.size() };
  
  typedef std::vector <
  std::shared_ptr<ComPWA::HelicityFormalism::SequentialTwoBodyDecay>::iterator seqDecayItr;
  
  seqDecayItr begin() { return _seqDecays.begin(); }
  
  seqDecayItr end() { return _seqDecays.end(); }
  
protected:
  virtual const double Integral();

  std::vector<std::shared_ptr<SequentialTwoBodyDecay> > _seqDecays;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_ */

/
