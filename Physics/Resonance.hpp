/*
 * Resonance.hpp
 *
 *  Created on: Mar 3, 2016
 *      Author: weidenka
 */

#ifndef CORE_RESONANCE_HPP_
#define CORE_RESONANCE_HPP_

#include <vector>
#include <memory>

#include <boost/iterator/filter_iterator.hpp>

#include "Core/AbsParameter.hpp"
#include "Core/Parameter.hpp"
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"

namespace ComPWA {

enum normStyle {
  none, /*!< no normaliztion between Amplitudes. */
  /*!< all amplitudes are normalized to one.
   *  The normalization factor is \f$ 1/\sqrt(\int |A|^2)\f$ */
  one
};

class Resonance {
public:
  Resonance(): _preFactor(1, 0){};
  virtual ~Resonance(){};

  //! Clone function
  virtual Resonance *Clone(std::string newName = "") const = 0;

  //! Get resonance name
  virtual std::string GetName() const { return _name; }

  //! Set resonance name
  virtual void SetName(std::string name) { _name = name; }

  //! Implementation of interface for streaming info about the strategy
  virtual std::string to_str() const = 0;

  //! Set prefactor
  virtual void SetPrefactor(std::complex<double> pre) { _preFactor = pre; }

  //! Get prefactor
  virtual std::complex<double> GetPrefactor() const { return _preFactor; }

  //! Get coefficient
  virtual std::complex<double> GetCoefficient() const {
    return std::polar(_magnitude->GetValue(), _phase->GetValue());
  }

  /**
   Get strength parameter

   @return strength parameter
   */
  std::shared_ptr<ComPWA::DoubleParameter> GetMagnitudePar() {
    return _magnitude;
  }

  /**
   Get strength parameter

   @return strength parameter
   */
  double GetMagnitude() const { return _magnitude->GetValue(); }

  /**
   Set strength parameter

   @param par Strength parameter
   */
  void SetMagnitudePar(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _magnitude = par;
  }

  /**
   Set strength parameter

   @param par Strength parameter
   */
  void SetMagnitude(double par) { _magnitude->SetValue(par); }

  /**
   Get phase parameter

   @return Phase parameter
   */
  std::shared_ptr<ComPWA::DoubleParameter> GetPhasePar() { return _phase; }

  /**
   Get phase parameter

   @return Phase parameter
   */
  double GetPhase() const { return _phase->GetValue(); }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  void SetPhasePar(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _phase = par;
  }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  void SetPhase(double par) { _phase->SetValue(par); }
  
  virtual void GetParameters(ParameterList& list) {
    list.AddParameter(GetMagnitudePar());
    list.AddParameter(GetPhasePar());
  }
  
  //! value of resonance at \param point
  virtual std::complex<double> Evaluate(const dataPoint &point) const = 0;

  //! Check of tree is available
  virtual bool HasTree() const { return false; }
  
  virtual std::shared_ptr<FunctionTree> GetTree(ParameterList &sample,
                                                  ParameterList &phspSample,
                                                  ParameterList &toySample,
                                                std::string suffix){
    return std::shared_ptr<FunctionTree>();
  }
  
  
protected:
  std::string _name;
  std::shared_ptr<ComPWA::DoubleParameter> _magnitude;
  std::shared_ptr<ComPWA::DoubleParameter> _phase;
  std::complex<double> _preFactor;
};
  
} /* namespace ComPWA */
#endif /* CORE_RESONANCE_HPP_ */
