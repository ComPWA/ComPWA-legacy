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
#include "Physics/AmplitudeSum/AmpWigner.hpp"

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
  void addBW(std::shared_ptr<AmpAbsDynamicalFunction>, std::shared_ptr<DoubleParameter>, std::shared_ptr<DoubleParameter>, std::shared_ptr<AmpWigner>);

  double evaluateSlice(std::complex<double>*, unsigned int, unsigned int) const ;
  double evaluate() const ;
  
  //double getMax() { return maxVal; };

protected:
  std::vector<std::shared_ptr<AmpAbsDynamicalFunction> > _pdfList ;   //  List of component PDFs
  std::vector<std::shared_ptr<DoubleParameter> > _intList;    //  List of relative intensities
  std::vector<std::shared_ptr<DoubleParameter> > _phaseList;  //  List of relative phases
  std::vector<std::shared_ptr<AmpWigner> > _angList ;   //  List of component angular distributions

  double maxVal;


};
 
#endif
