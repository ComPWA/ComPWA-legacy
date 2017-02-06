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
#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "Estimator/ChiOneD/ChiOneD.hpp"
#include "Core/ParameterList.hpp"

namespace ComPWA {

using DataReader::Data;
using ComPWA::ControlParameter;

namespace Estimator {
namespace ChiOneD {

ChiOneD::ChiOneD(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF) : pPIF_(inPIF), pDIF_(inDIF){

}

ChiOneD::~ChiOneD(){

}

std::shared_ptr<ControlParameter> ChiOneD::createInstance(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF){
  if(!instance_)
    instance_ = std::shared_ptr<ControlParameter>(new ChiOneD(inPIF, inDIF));

  return instance_;
}

double ChiOneD::controlParameter(ParameterList& minPar){
  unsigned int nBins = pDIF_->getNBins();

  double chi=0;
  for(unsigned int bin = 0; bin < nBins; bin++){
    double m12, weight;
    pDIF_->getBin(bin, m12, weight);

    std::vector<double> x;
    x.push_back(m12);
    ParameterList intensL = pPIF_->intensity(x);
    double intens = intensL.GetDoubleParameter(0)->GetValue();
    //double intens = pPIF_->intensity(x, minPar);

    chi += (weight - intens)*(weight - intens)/(double)nBins; //Just for binned data
  }

  return chi;
}

} /* namespace ChiOneD */
} /* namespace Estimator */
} /* namespace ComPWA */
