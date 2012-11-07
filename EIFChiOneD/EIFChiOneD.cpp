#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "EIFChiOneD.hpp"
#include "PWAEvent.hpp"
#include "PWAParticle.hpp"

EIFChiOneD::EIFChiOneD(std::shared_ptr<PIFBase> inPIF, std::shared_ptr<DIFBase> inDIF) : pPIF_(inPIF), pDIF_(inDIF){

}

EIFChiOneD::~EIFChiOneD(){

}

const int EIFChiOneD::getEstimatedVal(double& outChi){
  unsigned int nEvents = pDIF_->getNEvents();

  double chi=0;
  for(unsigned int evt = 0; evt < nEvents; evt++){
    PWAEvent theEvent;
    PWAParticle a, b;
    if(!(pDIF_->getEvent(evt, theEvent)))
      continue;
    theEvent.getParticle(0,a);
    theEvent.getParticle(1,b);
    double masssq = pow(a.getE()+b.getE(),2) - pow(a.getPx()+b.getPx() ,2) - pow(a.getPy()+b.getPy() ,2) - pow(a.getPz()+b.getPz() ,2);

    double intens = pPIF_->intensity(masssq, 1.5, 0.3); //TODO: par optimierbar

    chi += fabs(masssq - intens); //TODO: real Chi2
  }

  outChi = chi;

  return 0;
}
