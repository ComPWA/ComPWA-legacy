#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Core/PWAEvent.hpp"
#include "Core/PWAParticle.hpp"
#include "Core/ParameterList.hpp"

MinLogLH::MinLogLH(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF) : pPIF_(inPIF), pDIF_(inDIF){

}

MinLogLH::~MinLogLH(){

}

double MinLogLH::controlParameter(ParameterList& minPar){
  unsigned int nEvents = pDIF_->getNEvents();

  double norm = pPIF_->integral(minPar);

  double lh=0; //calculate LH:
  for(unsigned int evt = 0; evt < nEvents; evt++){
    PWAEvent theEvent;
    PWAParticle a, b;
    if( !(pDIF_->getEvent(evt, theEvent)) ){
      std::cout << "EIFChiOneD::controlParameter: Event not readable!" << std::endl; //TODO Exception
      continue;
    }
    if( !(theEvent.getParticle(0,a)) ){
      std::cout << "EIFChiOneD::controlParameter: Particle A not readable!" << std::endl; //TODO Exception
      continue;
    }
    if( !(theEvent.getParticle(1,b)) ){
      std::cout << "EIFChiOneD::controlParameter: Particle B not readable!" << std::endl; //TODO Exception
      continue;
    }

    double masssq = 0;
    masssq += pow(a.getE()+b.getE(),2);
    masssq -= pow(a.getPx()+b.getPx() ,2);
    masssq -= pow(a.getPy()+b.getPy() ,2);
    masssq -= pow(a.getPz()+b.getPz() ,2);

    std::vector<double> x;
    x.push_back(sqrt(masssq));
    double intens = pPIF_->intensity(x, minPar);

    lh -= (log(intens/norm));
  }

  return lh;
}
