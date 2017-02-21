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
  Resonance(){};
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
  virtual void SetPrefactor(std::complex<double> pre) = 0;

  //! Get prefactor
  virtual std::complex<double> GetPrefactor() const = 0;

  //! Get coefficient
  virtual std::complex<double> GetCoefficient() const = 0;

  //! Get magnitude of resonance name
  virtual double GetMagnitude() const = 0;

  //! Get magnitude of resonance id
  virtual std::shared_ptr<DoubleParameter> GetMagnitudePar() = 0;

  //! Get phase of resonance name
  virtual double GetPhase() const = 0;

  //! Get phase of resonance id
  virtual std::shared_ptr<DoubleParameter> GetPhasePar() = 0;

  //! value of resonance at \param point
  virtual std::complex<double> Evaluate(const dataPoint &point) const = 0;

  virtual std::shared_ptr<FunctionTree> SetupTree(ParameterList &sample,
                                                  ParameterList &toySample,
                                                  std::string suffix) = 0;
protected:
  std::string _name;
};
  
} /* namespace ComPWA */
#endif /* CORE_RESONANCE_HPP_ */
