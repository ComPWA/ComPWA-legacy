#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "Estimator/ChiOneD/ChiOneD.hpp"
#include "Core/ParameterList.hpp"

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
    ParameterList intensL = pPIF_->intensity(x, minPar);
    double intens = intensL.GetDoubleParameter(0)->GetValue();
    //double intens = pPIF_->intensity(x, minPar);

    chi += (weight - intens)*(weight - intens)/(double)nBins; //Just for binned data
  }

  return chi;
}
