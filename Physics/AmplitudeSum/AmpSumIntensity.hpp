//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     	Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding flatte type resonance, removing root
//dependence
//-------------------------------------------------------------------------------

#ifndef _AMPSUMINTENSITY_HPP
#define _AMPSUMINTENSITY_HPP

#include <vector>
#include <memory>
#include <map>
#include <string>

#include "Core/AmpIntensity.hpp"
#include "Core/Resonance.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Efficiency.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Generator.hpp"

#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/DalitzKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

  
class AmpSumIntensity : public AmpIntensity {

public:
  AmpSumIntensity(std::string name = "", normStyle ns = normStyle::none,
                  std::shared_ptr<Efficiency> eff =
                      std::shared_ptr<Efficiency>(new UnitEfficiency),
                  unsigned int nCalls = 30000);

  //! Copy constructor
  AmpSumIntensity(const AmpSumIntensity &copy);

  //! Destructor
  virtual ~AmpSumIntensity(){/* nothing */};

  //! Clone function
  virtual AmpSumIntensity *Clone(std::string newName = "") const;

  //===================== OPERATORS ==========================
  /** Operator for coherent addition of amplitudes
   *
   * @param other
   * @return
   */
  const AmpSumIntensity operator+(const AmpSumIntensity &other) const;

  /**! Add amplitude parameters to list
   * Add parameters only to list if not already in
   * @param list Parameter list to be filled
   */
  virtual void FillParameterList(ParameterList &list) const { };
  /** Operator for coherent addition of amplitudes
   *
   * @param rhs
   * @return
   */
  AmpSumIntensity &operator+=(const AmpSumIntensity &rhs);

  //===================== LOAD/SAVE CONFIG ====================
  //! Set prefactor
  virtual void SetPrefactor(std::complex<double> pre);

  //! Get Monte-Carlo precision (number of calls to integration algorithms)
  unsigned int GetMcPrecision() { return _nCalls; }

  //! get maximum value of amplitude with current parameters
  virtual double GetMaximum(std::shared_ptr<Generator> gen);

  //! Average width of all resonances
  virtual double averageWidth();

  //==================== PRINTING ============================
  //! Print overview over all amplitudes
  virtual void to_str();

  //! print all fit fractions; fitting errors are not available here
  virtual void printFractions();

  //==================== NORMALIZATION/INTEGRATION ==============
  //! Get normalization factor
  virtual double GetNormalization() const;

  //! Calculate inteference integral of two subresonances of @param A and @param
  //! B
  virtual double GetIntegralInterference(ampItr A, ampItr B) const;

  //! calculate integral for a list of resonances
  static double GetIntegralInterference(std::vector<ampItr> resList,
                                              unsigned int nCalls);

  /** Calculation of integral
   *
   * @param resList List of resonances to integrate
   * @param eff Efficiency correction (disable -> NULLPTR)
   * @param nCalls Monte-Carlo precision (number of calls)
   * @return
   */
  static double integral(std::vector<ampItr> resList,
                         std::shared_ptr<Efficiency> eff, int nCalls = 30000);

  //==================== EVALUATION =======================
  /**! Evaluate total amplitude
   * Using current set of parameters at phsp point @param point
   * @param point Point in phase-space
   * @return Complex function value
   */
  virtual std::complex<double> Evaluate(const dataPoint &point) const;

  /**! Evaluate total amplitude
   * Using current set of parameters at phsp point @param point . Amplitude is
   * multiplied with efficiency of datapoint.
   */
  virtual double Intensity(const dataPoint &point) const;

  /**! Evaluate total amplitude
   * Using current set of parameters at phsp point @param point .
   * No efficiency correction.
   */
  virtual double IntensityNoEff(const dataPoint &point) const;

  virtual double sliceIntensity(dataPoint &dataP, ParameterList &par,
                                      std::complex<double> *reso,
                                      unsigned int nResos, double N,
                                      unsigned int nF0, unsigned int nF2);

  //==================== FIT FRACTIONS =======================
  virtual void GetFitFractions(ParameterList &parList);

  static void GetFitFractions(ParameterList &parList, const AmpSumIntensity *amp);

  //================== ACCESS to resonances ====================
  //! Get ID of resonance from name
  virtual int GetIdOfResonance(std::string name);

  //! Get resonance name from ID
  virtual std::string GetNameOfResonance(unsigned int id);

  //! get resonance by @param name
  virtual std::shared_ptr<Amplitude> GetAmplitude(std::string name);

  //! get resonance by @param id
  virtual std::shared_ptr<Amplitude> GetAmplitude(unsigned int id);

  //! List of resonances (enabled AND disabled)
  virtual std::vector<std::shared_ptr<Amplitude>> &GetAmplitudes() {
    return _ampList;
  }

  virtual ampItr begin() {
    return _ampList.begin();
  }
  
  virtual ampItr end() {
    return _ampList.end();
  }

  //========== FunctionTree =============
  //! Check of tree is available
  virtual bool hasTree() { return 1; }
  //! Getter function for function tree
  virtual std::shared_ptr<FunctionTree> GetTree(ParameterList &sample,
                                                ParameterList &phspSample,
                                                ParameterList &toySample);

  //! Get amplitude integral for a set of (sub-) resonances
  virtual double Integral(std::vector<ampItr> resoList) const;
  
protected:
  //! Get amplitdude integral. Result is not efficiency corrected.
  virtual double Integral() const;
  
  //! Maximum value of amplitude. Necessary for event generation.
  double _maxFcnVal;
  //! Is amplitude maximum already calculated?
  bool _calcMaxFcnVal;
  //! calculate maximum value of amplitude with current parameters
  virtual void calcMaxVal(std::shared_ptr<Generator> gen);
  //! List of resonances
  std::vector<std::shared_ptr<Amplitude>> _ampList;
  //! Type of normalization
  normStyle _normStyle;
  //! precision for numeric integration
  unsigned int _nCalls;

  //========== FunctionTree =============
  /**Setup Basic Tree
   * @param sample sample of events for the amplitude calculation
   * @param phspSample sample of flat toy MC events for normalization of the
   * resonances
   * with efficiency corrected toy phsp sample or "normAcc" normalization tree
   * with sample
   * of accepted flat phsp events
   */
  std::shared_ptr<FunctionTree> setupBasicTree(ParameterList &sample,
                                               ParameterList &phspSample,
                                               std::string suffix = "");
  
};

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
