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

class SequentialTwoBodyDecay : public Amplitude {

public:
  //============ CONSTRUCTION ==================

  //! Function to create a full copy of the amplitude
  virtual Amplitude *Clone(std::string newName = "") const {
    auto tmp = (new SequentialTwoBodyDecay(*this));
    tmp->SetName(newName);
    return tmp;
  };
  
  /**
   Factory for SequentialTwoBodyDecay

   @param pt Configuration tree
   @return Constructed object
   */
  static std::shared_ptr<ComPWA::Physics::Amplitude>
  Factory(std::shared_ptr<Kinematics> kin, const boost::property_tree::ptree &pt);

  static boost::property_tree::ptree
  Save(std::shared_ptr<ComPWA::Physics::Amplitude> obj);

  //======= INTEGRATION/NORMALIZATION ===========
  
  //! Check of parameters have changed and normalization has to be recalculatecd
  bool CheckModified() const {
    if (Amplitude::CheckModified())
      return true;
    for (auto i : _partDecays)
      if (i->CheckModified())
        return true;
    return false;
  }
  
  //================ EVALUATION =================
  
  /**! Evaluate decay */
  virtual std::complex<double> Evaluate(const dataPoint &point) const {
    std::complex<double> result = GetCoefficient()*GetPreFactor();
    for (auto i : _partDecays)
      result *= i->Evaluate(point);

    assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
    return result;
  };

  //============ SET/GET =================

  /**
   Add a partial decay to Sequential decay

   @param d Partial decay
   */
  void Add(std::shared_ptr<ComPWA::Physics::Resonance> d) {
    _partDecays.push_back(d);
  }

  std::shared_ptr<ComPWA::Physics::Resonance> GetDecay(int pos) {
    return _partDecays.at(pos);
  }

  std::vector<std::shared_ptr<ComPWA::Physics::Resonance>> &GetDecays() {
    return _partDecays;
  }

  //! Fill ParameterList with fit fractions
  virtual void GetFitFractions(ParameterList &parList){};

  //! Fill ParameterList with fit fractions
  virtual void GetParameters(ParameterList &list);
  
  //! Fill vector with parameters
  virtual void GetParametersFast(std::vector<double> &list) const {
    Amplitude::GetParametersFast(list);
    for( auto i:_partDecays)
      i->GetParametersFast(list);
  }
  
  /*! Set phase space sample
   * We use the phase space sample to calculate the normalization. The sample
   * should be without efficiency applied.
   */
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample) {
    for (auto i : _partDecays)
      i->SetPhspSample(phspSample);
  }

  //======== ITERATORS/OPERATORS =============
  
  typedef std::vector<std::shared_ptr<ComPWA::Physics::Resonance>>::iterator
      partDecayItr;

  partDecayItr begin() { return _partDecays.begin(); }

  partDecayItr end() { return _partDecays.end(); }

  //=========== FUNCTIONTREE =================
  
  //! Check of tree is available
  virtual bool HasTree() const { return true; }

  /**! Setup function tree */
  virtual std::shared_ptr<FunctionTree> GetTree(std::shared_ptr<Kinematics> kin,
                                                const ParameterList &sample,
                                                const ParameterList &toySample,
                                                std::string suffix);

protected:
  std::vector<std::shared_ptr<ComPWA::Physics::Resonance>> _partDecays;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* SequentialTwoBodyDecay_h */
