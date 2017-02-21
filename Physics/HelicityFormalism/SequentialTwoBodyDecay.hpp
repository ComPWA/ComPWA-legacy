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

#include "Core/Amplitude.hpp"
#include "Core/Parameter.hpp"
#include "Physics/HelicityFormalism/PartialDecay.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class SequentialTwoBodyDecay : Amplitude {
  
public:
  /**! Evaluate decay */
  virtual std::complex<double> Evaluate(const dataPoint &point) const {
    std::complex<double> result =
        std::polar(_strength->GetValue(), _phase->GetValue());
    for( auto i : _partDecays)
      result *= i->Evaluate(point);
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
  Factory(const boost::property_tree::ptree &pt) {
    auto obj = std::shared_ptr<SequentialTwoBodyDecay>(); //std::make_shared(new SequentialTwoBodyDecay());
    obj->SetName( pt.get<string>("Amplitude.<xmlattr>.name","empty") );
    obj->SetStrength( ComPWA::DoubleParameterFactory(pt.get_child("strength")) );
    obj->SetPhase( ComPWA::DoubleParameterFactory(pt.get_child("phase")) );
    
    for(const auto& v : pt.get_child("Resonance") ){
      obj->Add( PartialDecay::Factory(v.second) );
    }
    return obj;
  }

  /**
   Add a partial decay to Sequential decay

   @param d Partial decay
   */
  void Add( std::shared_ptr<ComPWA::Physics::HelicityFormalism::PartialDecay> d ) {
    _partDecays.push_back(d);
  }
  
  /**
   Get strength parameter
   
   @return strength parameter
   */
  std::shared_ptr<ComPWA::DoubleParameter> GetStrength() {
    return _strength;
  }
  
  /**
   Get strength parameter
   
   @return strength parameter
   */
  double GetStrengthValue() { return _strength->GetValue(); }
  
  /**
   Set strength parameter
   
   @param par Strength parameter
   */
  void SetStrength(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _strength = par;
  }
  
  /**
   Set strength parameter
   
   @param par Strength parameter
   */
  void SetStrength(double par) { _strength->SetValue(par); }
  
  /**
   Get phase parameter
   
   @return Phase parameter
   */
  std::shared_ptr<ComPWA::DoubleParameter> GetPhase() {
    return _phase;
  }
  
  /**
   Get phase parameter
   
   @return Phase parameter
   */
  double GetPhaseValue() { return _phase->GetValue(); }
  
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

  
  /**
   Get number of partial decays

   @return Number of partial decays
   */
  size_t size() { return _partDecays.size(); };
  
  typedef std::vector<std::shared_ptr<PartialDecay> >::iterator partDecayItr;
  
  partDecayItr begin() { return _partDecays.begin(); }
  
  partDecayItr end() { return _partDecays.end(); }

protected:
  std::string name; // for full coherent amplitude construction reasons

  // TODO: we add this particle state info for the coherent sum stuff
  // the whole design is fucked up because of that, change that someday
  std::pair<ParticleStateInfo, std::pair<ParticleStateInfo, ParticleStateInfo>>
      decay_spin_info_;

  std::shared_ptr<DoubleParameter> _strength;
  std::shared_ptr<DoubleParameter> _phase;
  std::vector<std::shared_ptr<PartialDecay> > _partDecays;

};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* SequentialTwoBodyDecay_h */
