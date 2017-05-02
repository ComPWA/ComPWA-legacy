//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding UnitAmp
//-------------------------------------------------------------------------------
//! Physics Interface Base-Class.
/*! \class Amplitude
 * @file Amplitude.hpp
 * This class provides the interface to the model which tries to describe the
 * intensity. As it is pure virtual, one needs at least one implementation to
 * provide an model for the analysis which calculates intensities for an event
 * on
 * basis model parameters. If a new physics-model is derived from and fulfills
 * this base-class, no change in other modules are necessary to work with the
 * new
 * physics module.
 */

#ifndef AMPLITUDE_HPP_
#define AMPLITUDE_HPP_

#include <vector>
#include <memory>
#include <math.h>

//#include "Core/Resonance.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "DataReader/Data.hpp"

#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Generator.hpp"

namespace ComPWA {
namespace Physics {

class Amplitude {

public:
  //============ CONSTRUCTION ==================
  
  //! Constructor with an optional, unique name and an optional efficiency
  Amplitude(std::string name = "") : _name(name), _preFactor(1, 0) {}

  //! Destructor
  virtual ~Amplitude() { /* nothing */
  }

  //! Function to create a full copy of the amplitude
  virtual Amplitude *Clone(std::string newName = "") const = 0;

  //======= INTEGRATION/NORMALIZATION ===========
  
  //! Check of parameters have changed and normalization has to be recalculatecd
  bool CheckModified() const {
    if (GetMagnitude() != _current_magnitude || GetPhase() != _current_phase) {
      const_cast<double &>(_current_magnitude) = GetMagnitude();
      const_cast<double &>(_current_phase) = GetPhase();
      return true;
    }
    return false;
  }
  
  //================ EVALUATION =================
  
  /** Calculate value of amplitude at point in phase space
   *
   * @param point Data point
   * @return
   */
  virtual std::complex<double> Evaluate(const dataPoint &point) const = 0;

  //============ SET/GET =================
  
  //! Get name of amplitude
  virtual std::string GetName() const { return _name; }

  //! Set name of amplitude
  virtual void SetName(std::string name) { _name = name; }

  //! Get coefficient
  virtual std::complex<double> GetCoefficient() const {
    return std::polar(GetMagnitude(), GetPhase());
  }
  
  /** Update parameters
   *
   * @param par New list of parameters
   */
  virtual void UpdateParameters(ParameterList &par) { /* TODO */
  }

  //! Add parameters to list
  virtual void GetParameters(ParameterList &list) {
    list.AddParameter(_magnitude);
    list.AddParameter(_phase);
  }

  //! Fill ParameterList with fit fractions
  virtual void GetFitFractions(ParameterList &parList) = 0;
  
  /**
   Get Magnitude parameter

   @return Magnitude parameter
   */
  virtual std::shared_ptr<ComPWA::DoubleParameter> GetMagnitudeParameter() {
    return _magnitude;
  }

  /**
   Get Magnitude parameter

   @return Magnitude parameter
   */
  virtual double GetMagnitude() const {
    return std::fabs(_magnitude->GetValue());
  }

  /**
   Set Magnitude parameter

   @param par Magnitude parameter
   */
  virtual void
  SetMagnitudeParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _magnitude = par;
  }

  /**
   Set Magnitude parameter

   @param par Magnitude parameter
   */
  virtual void SetMagnitude(double par) { _magnitude->SetValue(par); }

  /**
   Get phase parameter

   @return Phase parameter
   */
  virtual std::shared_ptr<ComPWA::DoubleParameter> GetPhaseParameter() {
    return _phase;
  }

  /**
   Get phase parameter

   @return Phase parameter
   */
  virtual double GetPhase() const { return _phase->GetValue(); }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  virtual void SetPhaseParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _phase = par;
  }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  virtual void SetPhase(double par) { _phase->SetValue(par); }

  /**
 Set pre-factor

 @param par Pre-factor
 */
  virtual void SetPreFactor(std::complex<double> pre) { _preFactor = pre; }

  /**
 Get pre-factor

 @return Pre-factor
 */
  virtual std::complex<double> GetPreFactor() const { return _preFactor; }

  /*! Set phase space sample
   * We use the phase space sample to calculate the normalization. The sample
   * should be without efficiency applied.
   */
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample) = 0;

  //=========== FUNCTIONTREE =================
  
  //! Check of tree is available
  virtual bool HasTree() const { return 0; }

  //! Getter function for basic amp tree
  virtual std::shared_ptr<FunctionTree> GetTree(const ParameterList &sample,
                                                const ParameterList &toySample,
                                                std::string suffix) = 0;


protected:
  std::string _name;

  std::complex<double> _preFactor;

  std::shared_ptr<DoubleParameter> _magnitude;

  std::shared_ptr<DoubleParameter> _phase;

private:
  double _current_magnitude;
  double _current_phase;
};

typedef std::vector<std::shared_ptr<Amplitude>>::iterator ampItr;

} /* namespace Physics */
} /* namespace ComPWA */
#endif
