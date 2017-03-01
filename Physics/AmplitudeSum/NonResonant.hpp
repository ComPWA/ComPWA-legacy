/*
 * NonResonant.hpp
 *
 *  Created on: Jan 13, 2015
 *      Author: weidenka
 */

#ifndef NONRESONANT_HPP_
#define NONRESONANT_HPP_

#include "Core/DataPoint.hpp"
#include "Core/Parameter.hpp"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

class NonResonant : public AmpAbsDynamicalFunction {

public:
  NonResonant(normStyle nS = normStyle::one, int calls = 30000)
      : AmpAbsDynamicalFunction(nS, calls){};

  NonResonant(const char *name, std::shared_ptr<DoubleParameter> mag,
              std::shared_ptr<DoubleParameter> phase, std::string mother,
              std::string particleA, std::string particleB, int nCalls = 30000,
              normStyle nS = normStyle::one);

  //! Clone function
  virtual NonResonant *Clone(std::string newName = "") const {
    auto tmp = (new NonResonant(*this));
    if (newName != "")
      tmp->SetName(newName);
    return tmp;
  }

  //! Get resonance width
  virtual double GetWidth() const { return 0; }
  //! value of dynamical amplitude at \param point
  virtual std::complex<double> EvaluateAmp(const dataPoint &point) const{
    return dynamicalFunction();
  }
  //! value of WignerD amplitude at \param point
  virtual double EvaluateWignerD(const dataPoint &point) const { return 1; };

  //! Calculation integral |dynamical amplitude|^2
  virtual double GetIntegral() const{
    return Kinematics::instance()->GetPhspVolume();
  }

  static std::complex<double> dynamicalFunction();

  virtual std::shared_ptr<FunctionTree> GetTree(ParameterList &sample,
                                                  ParameterList &phspSample,
                                                  ParameterList &toySample,
                                                  std::string suffix);
  
  /** Get maximum value of amplitude
   * Maximum is numerically calculated using a random number generator
   * @param gen Random number generator
   * @return
   */
  virtual double GetMaxVal(std::shared_ptr<Generator> gen) { return 1; }
};
}
} /* namespace Physics */
} /* namespace ComPWA */
#endif /* NONRESONANT_HPP_ */
