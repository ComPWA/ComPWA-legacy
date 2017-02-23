//
//  Header.h
//  COMPWA
//
//  Created by Peter Weidenkaff on 21/02/2017.
//
//

#include "Core/AmpIntensity.hpp"
#include "Physics/HelicityFormalism/CoherentIntensity.hpp"

#ifndef INCOHERENT_INTENSITY_HPP
#define INCOHERENT_INTENSITY_HPP

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class IncoherentIntensity : public ComPWA::AmpIntensity {
public:
  //! Function to create a full copy of the amplitude
  AmpIntensity *Clone(std::string newName = "") const{
    auto tmp = (new IncoherentIntensity(*this));
    tmp->SetName(newName);
    return tmp;
  }
  
  /** Get maximum value of amplitude
   * Maximum is numerically calculated using a random number generator
   * @param gen Random number generator
   * @return
   */
  virtual double GetMaximum(std::shared_ptr<Generator> gen) {
    double max=0;
    for( auto i : _intens ){
      double tmp = i->GetMaximum(gen);
    if( max < tmp ) max=tmp;
    }
    return max;
  }
  
  //=========== EVALUATION =================
  /** Calculate intensity of amplitude at point in phase-space
   *
   * @param point Data point
   * @return
   */
  virtual double Intensity(const dataPoint &point) const {
    return ( IntensityNoEff(point)*_eff->evaluate(point) );
  };
  
  /** Calculate intensity of amplitude at point in phase-space
   * Intensity is calculated excluding efficiency correction
   * @param point Data point
   * @return
   */
  virtual double IntensityNoEff(const dataPoint &point) const {
    double result = 0;
    for( auto i : _intens ){
      result += i->IntensityNoEff(point);
    }
    return result;
  }
  
  //=========== PARAMETERS =================
  /**! Add amplitude parameters to list
   * Add parameters only to list if not already in
   * @param list Parameter list to be filled
   */
  virtual void FillParameterList(ParameterList &list) const { };
  
  //! Calculate & fill fit fractions of this amplitude to ParameterList
  virtual void GetFitFractions(ParameterList &parList) { };
  
  //! Getter function for basic amp tree
  virtual std::shared_ptr<FunctionTree>
  SetupTree(ParameterList &, ParameterList &, ParameterList &) {
    return std::shared_ptr<FunctionTree>();
  }
  
  /**
   Get number of partial decays
   
   @return Number of partial decays
   */
  size_t size(){ return _intens.size(); }
  
  typedef std::vector <
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::CoherentIntensity> >::
  iterator coherentIntItr;
  
  coherentIntItr begin() { return _intens.begin(); }
  
  coherentIntItr end() { return _intens.end(); }
  
  //============= PRINTING =====================
  //! Print amplitude to logging system
  virtual void to_str() { LOG(info) << "IncoherentIntensity"; }
  
protected:
  /** Calculate integral of amplitude.
   * The integral does not include efficiency correction
   */
  virtual const double Integral() {
    double result = 0;
    for( auto i : _intens ){
      /* Have to call GetNormalization() here since Integral() is protected */
      result += 1/i->GetNormalization();
    }
    return result;
  };
  
  std::vector<std::shared_ptr<ComPWA::Physics::HelicityFormalism::CoherentIntensity>> _intens;
  
};
  
} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
#endif /* Header_h */
