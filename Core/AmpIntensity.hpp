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
/*! \class AmpIntensity
 * @file AmpIntensity.hpp
 * This class provides the interface to the model which tries to describe the
 * intensity. As it is pure virtual, one needs at least one implementation to
 * provide an model for the analysis which calculates intensities for an event
 * on
 * basis model parameters. If a new physics-model is derived from and fulfills
 * this base-class, no change in other modules are necessary to work with the
 * new
 * physics module.
 */

#ifndef AMPINTENSITY_HPP_
#define AMPINTENSITY_HPP_

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
#include "Core/Amplitude.hpp"

namespace ComPWA {
  
class AmpIntensity {

public:
  //! Constructor with an optional, unique name and an optional efficiency
  AmpIntensity(std::string name = "",
            std::shared_ptr<Efficiency> eff =
                std::shared_ptr<Efficiency>(new UnitEfficiency))
      : _name(name), _eff(eff) {}

  //! Destructor
  virtual ~AmpIntensity() { /* nothing */
  }

  //! Function to create a full copy of the amplitude
  virtual AmpIntensity *Clone(std::string newName = "") const = 0;

  //============ SET/GET =================
  //! Get name of amplitude
  virtual std::string GetName() const { return _name; }

  //! Set name of amplitude
  virtual void SetName(std::string name) { _name = name; }

  //! Get efficiency
  virtual std::shared_ptr<Efficiency> GetEfficiency() { return _eff; };

  //! Set efficiency
  virtual void SetEfficiency(std::shared_ptr<Efficiency> eff) { _eff = eff; };

  /** Get maximum value of amplitude
   * Maximum is numerically calculated using a random number generator
   * @param gen Random number generator
   * @return
   */
  virtual double GetMaximum(std::shared_ptr<Generator> gen) const = 0;

  //============= PRINTING =====================
  //! Print amplitude to logging system
  virtual void to_str() const = 0;

  //=========== INTEGRATION/NORMALIZATION =================
  /** Calculate normalization of amplitude.
   * The integral includes efficiency correction
   */
  virtual double GetNormalization() const { return 1/Integral(); }

  //=========== EVALUATION =================
  /** Calculate intensity of amplitude at point in phase-space
   *
   * @param point Data point
   * @return
   */
  virtual double Intensity(const dataPoint &point) const = 0;

  /** Calculate intensity of amplitude at point in phase-space
   * Intensity is calculated excluding efficiency correction
   * @param point Data point
   * @return
   */
  virtual double IntensityNoEff(const dataPoint &point) const = 0;

  //! Add parameters to list
  virtual void GetParameters(ParameterList &list) const = 0;

  //! Fill ParameterList with fit fractions
  virtual void GetFitFractions(ParameterList &parList) = 0;

  //========== FUNCTIONTREE =============
  //! Check of tree is available
  virtual bool HasTree() { return 0; }

  //! Getter function for basic amp tree
  virtual std::shared_ptr<FunctionTree>
  GetTree(ParameterList &, ParameterList &, ParameterList &, std::string suffix="") {
    return std::shared_ptr<FunctionTree>();
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
  double GetStrengthValue() const { return _strength->GetValue(); }

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
  
protected:
  
  /** Calculate integral of amplitude.
   * The integral does not include efficiency correction
   */
  virtual double Integral() const = 0;

  //! Name
  std::string _name;

  //! Efficiency object
  std::shared_ptr<Efficiency> _eff;
  
  std::shared_ptr<ComPWA::DoubleParameter> _strength;
};
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
class GaussAmp : public AmpIntensity {
public:
  GaussAmp(const char *name, DoubleParameter _resMass,
           DoubleParameter _resWidth) {
    params.AddParameter(
        std::shared_ptr<DoubleParameter>(new DoubleParameter(_resMass)));
    params.AddParameter(
        std::shared_ptr<DoubleParameter>(new DoubleParameter(_resWidth)));
    initialise();
  }

  GaussAmp(const char *name, double _resMass, double _resWidth) {
    params.AddParameter(std::shared_ptr<DoubleParameter>(
        new DoubleParameter("mass", _resMass)));
    params.AddParameter(std::shared_ptr<DoubleParameter>(
        new DoubleParameter("width", _resWidth)));
    initialise();
  }

  //! Clone function
  GaussAmp *Clone(std::string newName = "") const {
    auto tmp = (new GaussAmp(*this));
    tmp->SetName(newName);
    return tmp;
  }

  virtual void initialise() {
    if (Kinematics::instance()->GetNVars() != 1)
      throw std::runtime_error("GaussAmp::initialize() | "
                               "this amplitude is for two body decays only!");
  };
  //! Clone function
  virtual GaussAmp *Clone(std::string newName = "") {
    auto tmp = (new GaussAmp(*this));
    tmp->SetName(newName);
    return tmp;
  }

  virtual void to_str() const {};

  virtual const double GetNormalization() { return 1/Integral(); }

  virtual double GetMaximum(std::shared_ptr<Generator> gen) const {
    double mass = params.GetDoubleParameter(0)->GetValue();
    std::vector<double> m;
    m.push_back(mass * mass);
    dataPoint p(m);
    return Intensity(p);
  }

  virtual void GetParameters(ParameterList &list) const { };
  
  virtual double Intensity(const dataPoint &point) const {

    double mass = params.GetDoubleParameter(0)->GetValue();
    double width = params.GetDoubleParameter(1)->GetValue();
    double sqrtS = std::sqrt(point.getVal(0));

    std::complex<double> gaus(
        std::exp(-1 * (sqrtS - mass) * (sqrtS - mass) / width / width / 2.), 0);

    return std::norm(gaus);
  }

  virtual double IntensityNoEff(const dataPoint &point) const {
    return Intensity(point);
  }

  virtual void GetFitFractions(ParameterList &parList) {}
  
protected:
  //! Get integral
  virtual double Integral() const{
    return (params.GetDoubleParameter(1)->GetValue() * std::sqrt(2 * M_PI));
  }
  
  //! List of interal parameters
  ParameterList params;
};
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/**! UnitAmp
 *
 * Example implementation of AmpIntensity with the function value 1.0 at all
 * points in PHSP. It is used to test the likelihood normalization.
 */
class UnitAmp : public AmpIntensity {
public:
  UnitAmp() {
    _eff = std::shared_ptr<Efficiency>(new UnitEfficiency());
  }

  virtual ~UnitAmp() { /* nothing */
  }

  virtual UnitAmp *Clone(std::string newName = "") const {
    auto tmp = new UnitAmp(*this);
    tmp->SetName(newName);
    return tmp;
  }

  virtual void to_str() const {}

  virtual double GetMaximum(std::shared_ptr<Generator> gen) const { return 1; }


  virtual const double GetNormalization() {
    LOG(info) << "UnitAmp::normalization() | "
                 "normalization not implemented!";
    return 1;
  }

  virtual double Intensity(const dataPoint &point) const {
    return _eff->Evaluate(point);
  }

  virtual void GetParameters(ParameterList &list) const { };
  
  virtual double IntensityNoEff(const dataPoint &point) const {
    return 1.0;
  }
  virtual void GetFitFractions(ParameterList &parList) {}

  //========== FunctionTree =============
  //! Check of tree is available
  virtual bool HasTree() { return 1; }

  //! Getter function for basic amp tree
  //! Getter function for basic amp tree
  virtual std::shared_ptr<FunctionTree>
  GetTree(ParameterList& sample, ParameterList& toySample, ParameterList& sample3) {
       return setupBasicTree(sample, toySample, "");
  }


protected:
  /**Setup Basic Tree
   *
   * @param sample data sample
   * @param toySample sample of flat toy MC events for normalization of the
   * resonances
   * @param suffix Which tree should be created? "data" data Tree, "norm"
   * normalization tree
   * with efficiency corrected toy phsp sample or "normAcc" normalization tree
   * with sample
   * of accepted flat phsp events
   */
  std::shared_ptr<FunctionTree> setupBasicTree(ParameterList &sample,
                                               ParameterList &toySample,
                                               std::string suffix) {

    int sampleSize = sample.GetMultiDouble(0)->GetNValues();

    LOG(debug) << "UnitAmp::setupBasicTree() generating new tree!";
    if (sampleSize == 0) {
      LOG(error) << "UnitAmp::setupBasicTree() data sample empty!";
      return std::shared_ptr<FunctionTree>();
    }
    std::shared_ptr<FunctionTree> newTree(new FunctionTree());
    // std::shared_ptr<MultAll> mmultDStrat(new MultAll(ParType::MDOUBLE));

    std::vector<double> oneVec(sampleSize, 1.0);
    std::shared_ptr<AbsParameter> one(new MultiDouble("one", oneVec));
    newTree->createHead("AmpIntensity" + suffix, one);
    std::cout << newTree->head()->to_str(10) << std::endl;
    return newTree;
  }
  
  virtual double Integral() const{
    return Kinematics::instance()->GetPhspVolume();
  }
};

} /* namespace ComPWA */
#endif
