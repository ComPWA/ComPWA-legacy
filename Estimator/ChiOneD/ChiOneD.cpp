#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "Estimator/ChiOneD/ChiOneD.hpp"
#include "Core/PWAParameter.hpp"

ChiOneD::ChiOneD(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF) : pPIF_(inPIF), pDIF_(inDIF){

}

ChiOneD::~ChiOneD(){

}

double ChiOneD::controlParameter(std::vector<std::shared_ptr<PWAParameter> >& minPar){
  unsigned int nBins = pDIF_->getNBins();

  double chi=0;
  for(unsigned int bin = 0; bin < nBins; bin++){
    double m12, weight;
    pDIF_->getBin(bin, m12, weight);

    std::vector<double> x;
    x.push_back(m12);
    double intens = pPIF_->intensity(x, minPar);

    chi += (weight - intens)*(weight - intens)/(double)nBins; //Just for binned data
  }

  return chi;
}
