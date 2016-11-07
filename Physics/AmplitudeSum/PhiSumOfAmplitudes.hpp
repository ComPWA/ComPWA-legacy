//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef PHISUMOFAMPLITUDES
#define PHISUMOFAMPLITUDES

#include "Core/Parameter.hpp"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

#include <vector>
#include <memory>

class TIterator;

namespace ComPWA {

class DoubleParameter;

namespace Physics {
namespace AmplitudeSum {
 
class PhiSumOfAmplitudes {
public:
  PhiSumOfAmplitudes() {} ; 
  PhiSumOfAmplitudes(const char *name);
  PhiSumOfAmplitudes(const PhiSumOfAmplitudes& other, const char* name=0) ;
//  virtual TObject* clone(const char* newname) const { return new PhiSumOfAmplitudes(*this,newname); }
  inline virtual ~PhiSumOfAmplitudes() { };

  void addBW(std::shared_ptr<AmpAbsDynamicalFunction> theRes , std::shared_ptr<DoubleParameter> r, std::shared_ptr<DoubleParameter> phi);
  void addBW(std::shared_ptr<AmpAbsDynamicalFunction>, std::shared_ptr<DoubleParameter>, std::shared_ptr<DoubleParameter>, std::shared_ptr<AmpWigner2>);
  
  virtual std::string GetName(){ return _name; };
  virtual std::string GetTitle(){ return GetName(); };
protected:

  double evaluate() const ;

private:
  std::string _name;

  std::vector<std::shared_ptr<AmpAbsDynamicalFunction> > 	_pdfList ;   //  List of component PDFs
  std::vector<std::shared_ptr<DoubleParameter> > 			_intList;    //  List of relative intensities
  std::vector<std::shared_ptr<DoubleParameter> > 			_phaseList;  //  List of relative phases
  std::vector<std::shared_ptr<AmpWigner2> >					_angList ;   //  List of component angular distributions

//  TIterator* _pdfIter  ;
//  TIterator* _intIter  ;
//  TIterator* _phaseIter;
//  TIterator* _angIter;

  //ClassDef(PhiSumOfAmplitudes,1) // Your description goes here...
};

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */
 
#endif
