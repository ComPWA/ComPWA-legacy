 
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

class CoherentIntensity : public AmpIntensity {

public:
  CoherentIntensity();
  virtual ~CoherentIntensity();

  virtual double GetMaximum(std::shared_ptr<Generator> gen);

  virtual double Intensity(const dataPoint &point) const;

  virtual double IntensityNoEff(const dataPoint &point) const;

  //========== FUNCTIONTREE =============
  //! Check of tree is available
  virtual bool hasTree() { return 0; }

  //! Getter function for basic amp tree
  virtual std::shared_ptr<FunctionTree>
  SetupTree(ParameterList &, ParameterList &, ParameterList &) {
    return std::shared_ptr<FunctionTree>();
  }

  AmpIntensity *Clone(std::string newName = "") const{
    auto tmp = (new CoherentIntensity(*this));
    tmp->SetName(newName);
    return tmp;
  }

  static std::shared_ptr<CoherentIntensity>
  Factory(const boost::property_tree::ptree &pt) {
    auto obj = std::shared_ptr<CoherentIntensity>();
    obj->SetName(pt.get<string>("AmpIntensity.<xmlattr>.name", "empty"));
    
    for(const auto& v : pt.get_child("Amplitude") ){
      obj->Add( SequentialTwoBodyDecay::Factory(v.second) );
    }
    return obj;
  }

  void Add(std::shared_ptr<SequentialTwoBodyDecay> d) {
    _seqDecays.push_back(d);
  }

  /**
   Get number of partial decays

   @return Number of partial decays
   */
  size_t size(){ return _seqDecays.size(); }

  typedef std::vector <
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::SequentialTwoBodyDecay> >::
          iterator seqDecayItr;

  seqDecayItr begin() { return _seqDecays.begin(); }

  seqDecayItr end() { return _seqDecays.end(); }
  
  //=========== PARAMETERS =================
  /**! Add amplitude parameters to list
   * Add parameters only to list if not already in
   * @param list Parameter list to be filled
   */
  virtual void FillParameterList(ParameterList &list) const {};
  
  //! Calculate & fill fit fractions of this amplitude to ParameterList
  virtual void GetFitFractions(ParameterList &parList) {};

  //============= PRINTING =====================
  //! Print amplitude to logging system
  virtual void to_str() { LOG(info) << "CoherentIntensity"; }
  
protected:
  virtual const double Integral();

  std::vector<std::shared_ptr<SequentialTwoBodyDecay>> _seqDecays;
};

//#define BOOST_TEST_MODULE CoherentIntensity_test
//#include <boost/test/unit_test.hpp>
//
//BOOST_AUTO_TEST_CASE( CoherentIntensity_test ) {
//  ptree tr;
//  stringstream ss;
//  read_xml(ss, tr);
//
//  BOOST_CHECK(CoherentIntensity::Factory(tr));
//}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_ */
