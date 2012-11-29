#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "EIFMinLogLH.hpp"
#include "PWAEvent.hpp"
#include "PWAParticle.hpp"
#include "PWAParameter.hpp"

EIFMinLogLH::EIFMinLogLH(std::shared_ptr<PIFBase> inPIF, std::shared_ptr<DIFBase> inDIF) : pPIF_(inPIF), pDIF_(inDIF){

}

EIFMinLogLH::~EIFMinLogLH(){

}

double EIFMinLogLH::controlParameter(std::vector<std::shared_ptr<PWAParameter> >& minPar){
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
