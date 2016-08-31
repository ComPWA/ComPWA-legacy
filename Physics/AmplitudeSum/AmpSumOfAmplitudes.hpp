//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - removing root dependence
//-------------------------------------------------------------------------------

#ifndef AMPSUMOFAMPLITUDES
#define AMPSUMOFAMPLITUDES

#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Core/DataPoint.hpp"

#include <vector>
#include <memory>

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

class AmpSumOfAmplitudes{
public:
  AmpSumOfAmplitudes();
  AmpSumOfAmplitudes(const char *name);
  AmpSumOfAmplitudes(const AmpSumOfAmplitudes& other, const char* name=0) ;
  virtual ~AmpSumOfAmplitudes();

  void addBW(std::shared_ptr<AmpAbsDynamicalFunction>, std::shared_ptr<DoubleParameter>, std::shared_ptr<DoubleParameter>);
  void addBW(std::shared_ptr<AmpAbsDynamicalFunction>, std::shared_ptr<DoubleParameter>, std::shared_ptr<DoubleParameter>, std::shared_ptr<AmpWigner2>);

  double evaluateSlice(dataPoint& point,std::complex<double>*, unsigned int, unsigned int, double, unsigned int, unsigned int) const ;
  std::complex<double> getFirstBW(dataPoint& point) const ;
  std::complex<double> getFirstReso(dataPoint& point) const ;
  std::complex<double> getFirstAmp(dataPoint& point) const ;
  double evaluate(const dataPoint& point) const ;
  
  //! Get resonance by name
  virtual std::shared_ptr<AmpAbsDynamicalFunction> getResonance(std::string name){
	  return _pdfList.at(getAmpId(name));
  }
  //! Get resonance by ID
  virtual std::shared_ptr<AmpAbsDynamicalFunction> getResonance(unsigned int id){
	  return _pdfList.at(id);
  }
  //! Get number of resonances
  virtual unsigned int getNumberOfResonances() {return _pdfList.size();};
  //! Get name of amplitude id
  virtual std::string getAmpName(unsigned int id) {return _pdfList.at(id)->GetName();};
  //! Get ID of amplitude name
  virtual int getAmpId(std::string name);
  //! Get intensity of amplitude id
  virtual double getAmpMagnitude(unsigned int id) { return _intList.at(id)->GetValue(); }
  //! Get phase of amplitude id
  virtual double getAmpPhase(unsigned int id) { return _phaseList.at(id)->GetValue(); }
  //! Get intensity of amplitude name
  virtual double getAmpMagnitude(std::string name) { return _intList.at(getAmpId(name))->GetValue(); }
  //! Get phase of amplitude id
  virtual double getAmpPhase(std::string name) { return _phaseList.at(getAmpId(name))->GetValue(); }
  //! Get spin of amplitude name
  virtual double getSpin(std::string name){
	  return getSpin(getAmpId(name));
  }
  //! Get spin of amplitude ID
  virtual double getSpin(unsigned int id);
  //! Get integral of amplitude name
  virtual double getAmpIntegral(std::string name){
	  return getAmpIntegral(getAmpId(name));
  }
  //! Get integral of amplitude ID
  virtual double getAmpIntegral(unsigned int id);
  //! Get |r|^2 Int(Amp) by name
  virtual double getAmpStrength(std::string name){
	  return getAmpStrength(getAmpId(name));
  }
  //! Get |r|^2 Int(Amp) by ID
  virtual double getAmpStrength(unsigned int id);

protected:
  //! Vector with resonances
  std::vector<std::shared_ptr<AmpAbsDynamicalFunction> > _pdfList ;   //  List of component PDFs
  //! Vector with magnitudes of resonances
  std::vector<std::shared_ptr<DoubleParameter> > _intList;    //  List of relative intensities
  //! Vector with phases of resonances
  std::vector<std::shared_ptr<DoubleParameter> > _phaseList;  //  List of relative phases
};
 
} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
