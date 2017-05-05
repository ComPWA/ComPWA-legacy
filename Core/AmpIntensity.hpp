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
//		 Peter Weidenkaff - adding UnitAmp
//-------------------------------------------------------------------------------

#ifndef AMPINTENSITY_HPP_
#define AMPINTENSITY_HPP_

#include <vector>
#include <memory>
#include <math.h>

#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "DataReader/Data.hpp"

#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Generator.hpp"

namespace ComPWA {

/*! \class AmpIntensity
 * This class provides the interface an amplitude intensity. The intensity can
 * be moduled with a (double) strength parameter. The intensity has
 * to be normalized to one when integrated over the phase space.
 * Since the intensity is a physically observable quantity it has to be
 * corrected for the space depended reconstruction efficiency. The normalization
 * has to take this into account as well.
 */
class AmpIntensity {

public:
  //============ CONSTRUCTION ==================

  //! Constructor with an optional name, strength and efficiency
  AmpIntensity(std::string name = "",
               std::shared_ptr<DoubleParameter> strength =
                   std::shared_ptr<DoubleParameter>(new DoubleParameter("",
                                                                        1.0)),
               std::shared_ptr<Efficiency> eff =
                   std::shared_ptr<Efficiency>(new UnitEfficiency))
      : _name(name), _eff(eff), _strength(strength), _current_strength(-999) {
    _strength->FixParameter(true);
  }

  virtual ~AmpIntensity() { /* nothing */
  }

  //! Function to create a full copy
  virtual AmpIntensity *Clone(std::string newName = "") const = 0;

  //======= INTEGRATION/NORMALIZATION ===========

  //! Calculate normalization
  virtual double GetNormalization() const { return 1.0 / Integral(); }

  //! Check if parameters have changed
  bool CheckModified() const {
    if (GetStrength() != _current_strength) {
      const_cast<double &>(_current_strength) = GetStrength();
      return true;
    }
    return false;
  }

  //================ EVALUATION =================

  /*! Evaluate intensity at dataPoint in phase-space
   * @param point Data point
   * @return Intensity
   */
  virtual double Intensity(const dataPoint &point) const = 0;

  /*! Evaluate intensity at dataPoint in phase-space (excluding normalization).
   * @param point Data point
   * @return Intensity
   */
  virtual double IntensityNoNorm(const dataPoint &point) const {return 1.0;};

  //============ SET/GET =================
  //! Get name
  virtual std::string GetName() const { return _name; }

  //! Set name
  virtual void SetName(std::string name) { _name = name; }

  //! Get efficiency
  virtual std::shared_ptr<Efficiency> GetEfficiency() { return _eff; };

  //! Set efficiency
  virtual void SetEfficiency(std::shared_ptr<Efficiency> eff) { _eff = eff; };

  //! Get strength parameter
  std::shared_ptr<ComPWA::DoubleParameter> GetStrengthParameter() {
    return _strength;
  }

  //! Get strength parameter
  double GetStrength() const { return _strength->GetValue(); }

  //! Set strength parameter
  void SetStrengthParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _strength = par;
  }

  //! Set strength parameter
  void SetStrength(double par) { _strength->SetValue(par); }

  /*! Get maximum value of amplitude.
   * Maximum is numerically calculated using a random number generator
   */
  virtual double GetMaximum(std::shared_ptr<Generator> gen) const = 0;

  virtual void GetParameters(ParameterList &list) = 0;
  
  //! Fill vector with parameters
  virtual void GetParametersFast(std::vector<double> &list) const {
    list.push_back(GetStrength());
  }

  //! Fill ParameterList with fit fractions
  virtual void GetFitFractions(ParameterList &parList) = 0;

  /*! Set phase space samples
   * We use phase space samples to calculate the normalizations. In case of
   * intensities we phase space sample phspSample includes the event efficiency.
   * The sample toySample is used for normalization calculation for e.g.
   * Resonacnes without efficiency.
   */
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::dataPoint>> toySample) = 0;

  virtual std::shared_ptr<AmpIntensity> GetComponent(std::string name){};
  
  virtual void Reset() {};
  //========== FUNCTIONTREE =============

  //! Check of tree is available
  virtual bool HasTree() const { return false; }

  /*! Get FunctionTree
   * @param sample Data sample
   * @param phspSample Sample of phase space distributed events including
   * efficiency.
   * @param toySample Sample of phase space distributed events without
   * efficiency.
   */
  virtual std::shared_ptr<FunctionTree> GetTree(const ParameterList &sample,
                                                const ParameterList &phspSample,
                                                const ParameterList &toySample,
                                                std::string suffix = "") {
    return std::shared_ptr<FunctionTree>();
  }

  //======== ITERATORS/OPERATORS =============

protected:
  /*! Calculate integral.
   * Since AmpIntensity represents an intensity the integration has to
   * incorporate the phase space depended efficiency.
   */
  virtual double Integral() const = 0;

  //! Name
  std::string _name;

  //! Phase space depended efficiency
  std::shared_ptr<Efficiency> _eff;

  std::shared_ptr<ComPWA::DoubleParameter> _strength;

private:
  //! temporary strength
  double _current_strength;
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
    if (Kinematics::Instance()->GetNVars() != 1)
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

  virtual double GetNormalization() const { return 1 / Integral(); }

  virtual double GetMaximum(std::shared_ptr<Generator> gen) const {
    double mass = params.GetDoubleParameter(0)->GetValue();
    std::vector<double> m;
    m.push_back(mass * mass);
    dataPoint p(m);
    return Intensity(p);
  }

  virtual void GetParameters(ParameterList &list){};

  virtual double Intensity(const dataPoint &point) const {

    double mass = params.GetDoubleParameter(0)->GetValue();
    double width = params.GetDoubleParameter(1)->GetValue();
    double sqrtS = std::sqrt(point.GetValue(0));

    std::complex<double> gaus(
        std::exp(-1 * (sqrtS - mass) * (sqrtS - mass) / width / width / 2.), 0);

    return std::norm(gaus);
  }

  virtual double IntensityNoNorm(const dataPoint &point) const {
    return Intensity(point);
  }

  virtual void GetFitFractions(ParameterList &parList) {}

  /*! Set phase space samples
   * We use phase space samples to calculate the normalizations. In case of
   * intensities we phase space sample phspSample includes the event efficiency.
   * The sample toySample is used for normalization calculation for e.g.
   * Resonacnes without efficiency.
   */
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::dataPoint>> toySample) {}

protected:
  //! Get integral
  virtual double Integral() const {
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
  UnitAmp() { _eff = std::shared_ptr<Efficiency>(new UnitEfficiency()); }

  virtual ~UnitAmp() { /* nothing */
  }

  virtual UnitAmp *Clone(std::string newName = "") const {
    auto tmp = new UnitAmp(*this);
    tmp->SetName(newName);
    return tmp;
  }

  virtual void to_str() const {}

  virtual double GetMaximum(std::shared_ptr<Generator> gen) const { return 1; }

  virtual double GetNormalization() const {
    LOG(info) << "UnitAmp::normalization() | "
                 "normalization not implemented!";
    return 1;
  }

  virtual double Intensity(const dataPoint &point) const {
    return _eff->Evaluate(point);
  }

  virtual void GetParameters(ParameterList &list){};

  virtual double IntensityNoNorm(const dataPoint &point) const { return 1.0; }
  
  virtual void GetFitFractions(ParameterList &parList) {}

  //========== FunctionTree =============
  //! Check of tree is available
  virtual bool HasTree() const { return 1; }

  //! Getter function for basic amp tree
  //! Getter function for basic amp tree
  virtual std::shared_ptr<FunctionTree> GetTree(const ParameterList &sample,
                                                const ParameterList &toySample,
                                                const ParameterList &sample3,
                                                std::string suffix = "") {
    return setupBasicTree(sample, toySample, "");
  }

  /*! Set phase space samples
   * We use phase space samples to calculate the normalizations. In case of
   * intensities we phase space sample phspSample includes the event efficiency.
   * The sample toySample is used for normalization calculation for e.g.
   * Resonacnes without efficiency.
   */
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample,
                std::shared_ptr<std::vector<ComPWA::dataPoint>> toySample) {}

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
  std::shared_ptr<FunctionTree> setupBasicTree(const ParameterList &sample,
                                               const ParameterList &toySample,
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

  virtual double Integral() const {
    return Kinematics::Instance()->GetPhspVolume();
  }
};

  //! Split string into pieces which are separated by blanks
inline std::vector<std::string> splitString(std::string str) {
  std::vector<std::string> result;
  std::istringstream iStr(str);
  std::vector<std::string> stringFrag{std::istream_iterator<std::string>{iStr},
                                      std::istream_iterator<std::string>{}};
  for (auto i : stringFrag) {
    result.push_back(i);
  }
  return result;
}

} /* namespace ComPWA */
#endif
