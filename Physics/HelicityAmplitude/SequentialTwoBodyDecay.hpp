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

#include "Core/Amplitude.hpp"
#include "Physics/DynamicalDecayFunctions/AbstractDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class SequentialTwoBodyDecay : Amplitude {

  
public:
  /**! Evaluate decay */
  std::complex<double> Evaluate(dataPoint &point) {
    std::complex<double> result =
        std::polar(_strength->GetValue(), _phase->GetValue());
    for( auto i : _partDecays)
      result *= (*i)->Evaluate(point);
    return result;
  };

  /**! Setup function tree */
  virtual std::shared_ptr<FunctionTree> SetupTree(ParameterList &sample,
                                                  ParameterList &toySample,
                                                  std::string suffix) {
    return std::shared_ptr<FunctionTree>();
  };

  /**
   Factory for SequentialTwoBodyDecay

   @param pt Configuration tree
   @return Constructed object
   */
  static std::shared_ptr<SequentialTwoBodyDecay>
  Factory(boost::property_tree::ptree &pt) {
    auto obj = std::make_shared<SeqentialTwoBodyDecay>();
    obj->SetName( pt.get<string>("Amplitude.<xmlattr>.name","empty") );
    obj->SetStrength(DoubleParameterFactory(pt.get_child("strength")));
    obj->SetPhase(DoubleParameterFactory(pt.get_child("phase")));
    BOOST_FOREACH (ptree::value_type const &ampTree,
                   pt.second.get_child("Resonance")) {
      obj->Add( PartialDecay::Factory(ampTree) );
    }
    return obj;
  }

  /**
   Add a partial decay to Sequential decay

   @param d Partial decay
   */
  void
  Add( std::shared_ptr<ComPWA::HelicityFormalism::PartialDecay> d ) {
    _partDecays.push_back(d);
  }
  
  /**
   Get number of partial decays

   @return Number of partial decays
   */
  size_t size() { return _partDecays.size() };
  
  typedef std::vector <
  std::shared_ptr<ComPWA::HelicityFormalism::PartialDecay>::iterator partDecayItr;
  
  partDecayItr begin() { return _partDecays.begin(); }
  
  partDecayItr end() { return _partDecays.end(); }

protected:
  std::string name; // for full coherent amplitude construction reasons

  // TODO: we add this particle state info for the coherent sum stuff
  // the whole design is fucked up because of that, change that someday
  std::pair<ParticleStateInfo, std::pair<ParticleStateInfo, ParticleStateInfo>>
      decay_spin_info_;

  std::shared_ptr<ComPWA::DoubleParameter> _strength;
  std::shared_ptr<ComPWA::DoubleParameter> _phase;
  std::vector <
      std::shared_ptr<ComPWA::HelicityFormalism::PartialDecay> _partDecays;

  SequentialTwoBodyDecayAmp
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDE_HPP_ */

/

#endif /* SequentialTwoBodyDecay_h */
