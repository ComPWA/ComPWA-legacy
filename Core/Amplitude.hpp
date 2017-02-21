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

#include "Core/Resonance.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "DataReader/Data.hpp"

#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Generator.hpp"

namespace ComPWA {

class Amplitude {

public:
  //! Constructor with an optional, unique name and an optional efficiency
  Amplitude(std::string name = "",
            std::shared_ptr<Efficiency> eff =
                std::shared_ptr<Efficiency>(new UnitEfficiency))
      : _name(name), eff_(eff) {}

  //! Destructor
  virtual ~Amplitude() { /* nothing */
  }

  //! Function to create a full copy of the amplitude
  virtual Amplitude *Clone(std::string newName = "") const = 0;

  //============ SET/GET =================
  //! Get name of amplitude
  virtual std::string GetName() const { return _name; }

  //! Set name of amplitude
  virtual void SetName(std::string name) { _name = name; }

  //! Get efficiency
  virtual std::shared_ptr<Efficiency> GetEfficiency() { return eff_; };

  //! Set efficiency
  virtual void SetEfficiency(std::shared_ptr<Efficiency> eff) { eff_ = eff; };

  /** Get maximum value of amplitude
   * Maximum is numerically calculated using a random number generator
   * @param gen Random number generator
   * @return
   */
  virtual double GetMaxVal(std::shared_ptr<Generator> gen) = 0;
  //	virtual bool copyParameterList(ParameterList& par) =0;

  //============= PRINTING =====================
  //! Print amplitude to logging system
  virtual void to_str() = 0;

  //=========== INTEGRATION/NORMALIZATION =================
  /** Calculate normalization of amplitude.
   * The integral includes efficiency correction
   */
  virtual const double GetNormalization() = 0;

  /** Calculate integral of amplitude.
   * The integral does not include efficiency correction
   */
  virtual const double GetIntegral() = 0;

  /** Calculate integral of amplitudes for a given vector.
   * The integral does not include efficiency correction
   *
   * @param resoList Vector of amplitudes
   * @return
   */
//  virtual const double GetIntegral(std::vector<resonanceItr> resoList) {
//    return 0;
//  };

  /** Calculate interference integral
   *
   * @param A First resonance
   * @param B Second resonance
   * @return a*conj(b)+conj(a)*b
   */
//  virtual const double GetIntegralInterference(resonanceItr A, resonanceItr B) {
//    return -999;
//  };

  //=========== EVALUATION =================
  /** Calculate value of amplitude at point in phase space
   *
   * @param point Data point
   * @return
   */
 
  virtual std::complex<double> Evaluate(const dataPoint &point) const = 0;

  //=========== PARAMETERS =================
  /** Update parameters of resonance
   *
   * @param par New list of parameters
   */
  virtual void UpdateParameters(ParameterList &par);

  /**! Add amplitude parameters to list
   * Add parameters only to list if not already in
   * @param list Parameter list to be filled
   */
  virtual void FillParameterList(ParameterList &list) const;

  //! Calculate & fill fit fractions of this amplitude to ParameterList
  virtual void GetFitFractions(ParameterList &parList) = 0;

  //============= ACCESS TO RESONANCES ================
//  //! Iterator on first resonance (which is enabled)
//  virtual resonanceItr GetResonanceItrFirst(){};
//
//  //! Iterator on last resonance (which is enabled)
//  virtual resonanceItr GetResonanceItrLast(){};
//
//  //! Iterator on last resonance (which is enabled)
//  virtual const std::vector<resonanceItr> GetResonanceItrList() {
//    return std::vector<resonanceItr>();
//  };

  //========== FUNCTIONTREE =============
  //! Check of tree is available
  virtual bool hasTree() { return 0; }

  //! Getter function for basic amp tree
  virtual std::shared_ptr<FunctionTree>
  GetTree(ParameterList &, ParameterList &, ParameterList &) {
    return std::shared_ptr<FunctionTree>();
  }

  /**! Check if amplitudes have a FunctionTree
   *
   * @param ampVec corresponding vector of amplitudes
   * @return true if tree is present
   */
  static bool AmpHasTree(std::vector<std::shared_ptr<Amplitude>> ampVec);

  /**! Return Precision of Monte-Carlo Normalization
   *
   * @return integer value of used precision
   */
  unsigned int GetMcPrecision() { return 30000; } // Todo: fixed?

protected:
  //! Name
  std::string _name;

  //! need to store this object for boost::filter_iterator
//  resIsEnabled _resEnabled;

  //! Amplitude value
  ParameterList result;

  //! List of interal parameters
  ParameterList params;

  //! Efficiency object
  std::shared_ptr<Efficiency> eff_;
};

} /* namespace ComPWA */
#endif
