//
//  Header.h
//  COMPWA
//
//  Created by Peter Weidenkaff on 21/02/2017.
//
//

#include "Core/AmpIntensity.hpp"
#include "Physics/HelicityFormalism/CoherentIntensity.hpp"
#include "Tools/Integration.hpp"

#ifndef INCOHERENT_INTENSITY_HPP
#define INCOHERENT_INTENSITY_HPP

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

  

class IncoherentIntensity : public ComPWA::AmpIntensity {

public:
  //============ CONSTRUCTION ==================

  IncoherentIntensity() : ComPWA::AmpIntensity() {}

  //! Function to create a full copy of the amplitude
  ComPWA::AmpIntensity *Clone(std::string newName = "") const {
    auto tmp = (new IncoherentIntensity(*this));
    tmp->SetName(newName);
    return tmp;
  }

  static std::shared_ptr<IncoherentIntensity>
  Factory(const boost::property_tree::ptree &pt);

  static boost::property_tree::ptree
  Save(std::shared_ptr<IncoherentIntensity> intens);

  //================ EVALUATION =================

  /** Calculate intensity of amplitude at point in phase-space
   *
   * @param point Data point
   * @return
   */
  virtual double Intensity(const ComPWA::dataPoint &point) const {
    
    // We have to get around the constness of the interface definition.
    std::vector<std::vector<double>> parameters(_parameters);
    std::vector<double> normValues(_normValues);
    if(_intens.size() != parameters.size() )
      parameters = std::vector<std::vector<double>>(_intens.size());
    if(_intens.size() != normValues.size() )
      normValues = std::vector<double>(_intens.size());
    
    double result = 0;
    for(int i=0; i< _intens.size(); i++){
      std::vector<double> params;
      _intens.at(i)->GetParametersFast(params);
      if ( parameters.at(i) != params ) { //recalculate normalization
        parameters.at(i) = params;
        normValues.at(i) = 1/Tools::Integral(_intens.at(i), _phspSample);
        normValues.at(i) *= _intens.at(i)->GetStrength();
      }
      result += _intens.at(i)->Intensity(point) * normValues.at(i);
    }
    
    const_cast<std::vector<std::vector<double>>&>(_parameters) = parameters;
    const_cast<std::vector<double>&>(_normValues) = normValues;

    return (GetStrength() * result);
  }

  //================== SET/GET =================

  void AddIntensity(std::shared_ptr<ComPWA::AmpIntensity> intens) {
    _intens.push_back(intens);
  }

  std::shared_ptr<ComPWA::AmpIntensity> GetIntensity(int pos) {
    return _intens.at(pos);
  }

  std::vector<std::shared_ptr<ComPWA::AmpIntensity>> &GetIntensities() {
    return _intens;
  }

  virtual void Reset() {
    _intens.clear();
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
    for (auto i : _intens) {
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
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::dataPoint>> toySample) {
    _phspSample = phspSample;
    for (auto i : _intens) {
      i->SetPhspSample(phspSample, toySample);
    }
  };

  virtual std::shared_ptr<AmpIntensity> GetComponent(std::string name);

  //======== ITERATORS/OPERATORS =============
  typedef std::vector<std::shared_ptr<ComPWA::AmpIntensity>>::iterator
      coherentIntItr;

  coherentIntItr First() { return _intens.begin(); }

  coherentIntItr Last() { return _intens.end(); }

  //=========== FUNCTIONTREE =================

  //! Check of tree is available
  virtual bool HasTree() const { return true; }

  //! Get FunctionTree
  virtual std::shared_ptr<ComPWA::FunctionTree>
  GetTree(const ComPWA::ParameterList &sample,
          const ComPWA::ParameterList &phspSample,
          const ComPWA::ParameterList &toySample, std::string suffix = "");

protected:
  //! Phase space sample to calculate the normalization and maximum value.
  std::shared_ptr<std::vector<ComPWA::dataPoint>> _phspSample;
  
  // Caching of normalization values
  std::vector<double> _normValues;
  
  // Temporary storage of the para
  std::vector<std::vector<double>> _parameters;

  std::vector<std::shared_ptr<ComPWA::AmpIntensity>> _intens;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
#endif /* Header_h */
