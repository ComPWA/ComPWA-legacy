// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef GAUSSAMP_HPP_
#define GAUSSAMP_HPP_

#include <vector>
#include <memory>
#include <math.h>

#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "DataReader/Data.hpp"

#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Generator.hpp"
#include "Core/AmpIntensity.hpp"

namespace ComPWA {

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
    tmp->_name(newName);
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
    tmp->_name(newName);
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

  virtual std::shared_ptr<AmpIntensity> GetComponent(std::string name) {
    return std::shared_ptr<AmpIntensity>();
  }

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
    tmp->_name(newName);
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

  virtual std::shared_ptr<AmpIntensity> GetComponent(std::string name) {
    return std::shared_ptr<AmpIntensity>();
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

} /* namespace ComPWA */


#endif /* CORE_GAUSSAMP_HPP_ */
