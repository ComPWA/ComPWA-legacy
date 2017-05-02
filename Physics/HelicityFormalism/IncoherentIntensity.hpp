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

  //======= INTEGRATION/NORMALIZATION ===========

  //! Check if parameters of this class or one of its members have changed
  bool CheckModified() const {
    if (AmpIntensity::CheckModified())
      return true;
    for (auto i : _intens)
      if (i->CheckModified())
        return true;

    return false;
  }

  //================ EVALUATION =================

  /** Calculate intensity of amplitude at point in phase-space
   *
   * @param point Data point
   * @return
   */
  virtual double Intensity(const ComPWA::dataPoint &point) const {
    return (IntensityNoEff(point) * _eff->Evaluate(point));
  };

  /** Calculate intensity of amplitude at point in phase-space
   * Intensity is calculated excluding efficiency correction
   * @param point Data point
   * @return
   */
  virtual double IntensityNoEff(const ComPWA::dataPoint &point) const {
    double result = 0;
    for (auto i : _intens) {
      result += i->IntensityNoEff(point);
    }
    return GetStrength() * result;
  }

  //================== SET/GET =================

  /** Get maximum value of amplitude
   * We ask for the maximum of the coherent intensities and use the largest
   * value
   * @param gen Random number generator
   * @return
   */
  virtual double GetMaximum(std::shared_ptr<ComPWA::Generator> gen) const {
    double max = 0;
    for (auto i : _intens) {
      double tmp = i->GetMaximum(gen);
      if (max < tmp)
        max = tmp;
    }
    return max;
  }

  void
  Add(std::shared_ptr<ComPWA::Physics::HelicityFormalism::CoherentIntensity>
          d) {
    _intens.push_back(d);
  }

  std::shared_ptr<ComPWA::Physics::HelicityFormalism::CoherentIntensity>
  GetIntensity(int pos) {
    return _intens.at(pos);
  }

  std::vector<
      std::shared_ptr<ComPWA::Physics::HelicityFormalism::CoherentIntensity>> &
  GetIntensities() {
    return _intens;
  }

  /**! Add amplitude parameters to list
   * Add parameters only to list if not already in
   * @param list Parameter list to be filled
   */
  virtual void GetParameters(ComPWA::ParameterList &list);

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
    for (auto i : _intens) {
      i->SetPhspSample(phspSample, toySample);
    }
  };
  
  //======== ITERATORS/OPERATORS =============
  typedef std::vector<std::shared_ptr<
      ComPWA::Physics::HelicityFormalism::CoherentIntensity>>::iterator
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
  /*! Calculate integral of amplitude.
   * CoherentIntensities are required to be normalized. Therefore, we
   * sum up the strength's to get the integral.
   */
  virtual double Integral() const {
    double fractionSum = 0;
    for (auto i : _intens) {
      fractionSum += i->GetStrength();
    }
    return fractionSum;
  };

  std::vector<
      std::shared_ptr<ComPWA::Physics::HelicityFormalism::CoherentIntensity>>
      _intens;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
#endif /* Header_h */
