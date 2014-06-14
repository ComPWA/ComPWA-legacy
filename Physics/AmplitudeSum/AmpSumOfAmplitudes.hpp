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

//class TIterator;
//class RooRealVar;
 
class AmpSumOfAmplitudes{
public:
  AmpSumOfAmplitudes();
  AmpSumOfAmplitudes(const char *name);
  AmpSumOfAmplitudes(const AmpSumOfAmplitudes& other, const char* name=0) ;
  virtual ~AmpSumOfAmplitudes();

  void addBW(std::shared_ptr<AmpAbsDynamicalFunction>, std::shared_ptr<DoubleParameter>, std::shared_ptr<DoubleParameter>);
  void addBW(std::shared_ptr<AmpAbsDynamicalFunction>, std::shared_ptr<DoubleParameter>, std::shared_ptr<DoubleParameter>, std::shared_ptr<AmpWigner2>);

  double evaluateSlice(dataPoint& point,std::complex<double>*, unsigned int, unsigned int) const ;
  std::complex<double> getFirstBW(dataPoint& point) const ;
  std::complex<double> getFirstReso(dataPoint& point) const ;
  std::complex<double> getFirstAmp(dataPoint& point) const ;
  double evaluate(dataPoint& point) const ;
  
  virtual double getTotalIntegral(std::string name);
  virtual double getTotalIntegral(unsigned int id);
  virtual double getUnormalizedFraction(std::string name);
  virtual double getUnormalizedFraction(unsigned int id);
  virtual std::shared_ptr<AmpAbsDynamicalFunction> getResonance(std::string name);
  virtual std::shared_ptr<AmpAbsDynamicalFunction> getResonance(unsigned int id){ return _pdfList[id]; }
  virtual unsigned int getNAmps() {return _pdfList.size();};
  virtual std::string getAmpName(unsigned int id) {return _pdfList[id]->GetName();};
protected:
  std::vector<std::shared_ptr<AmpAbsDynamicalFunction> > _pdfList ;   //  List of component PDFs
  std::vector<std::shared_ptr<DoubleParameter> > _intList;    //  List of relative intensities
  std::vector<std::shared_ptr<DoubleParameter> > _phaseList;  //  List of relative phases
  std::vector<std::shared_ptr<AmpWigner2> > _angList ;   //  List of component angular distributions

  double maxVal;


};
 
#endif
