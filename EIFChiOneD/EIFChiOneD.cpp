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

double EIFChiOneD::controlParameter(const std::vector<double>& minPar){
  unsigned int nEvents = pDIF_->getNEvents();

  double chi=0;
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
    //const double ta(a.getE());
   // std::cout << "Test" << std::endl;
    //std::cout << "Event, E, Px: \t" << evt << "\t" << a.getE() << "\t" << a.getPx() << std::endl;
    double masssq = 0;
    masssq += pow(a.getE()+b.getE(),2);
    masssq -= pow(a.getPx()+b.getPx() ,2);
    masssq -= pow(a.getPy()+b.getPy() ,2);
    masssq -= pow(a.getPz()+b.getPz() ,2);

    //std::cout << "Event, Par0, Par1: \t" << evt << "\t" << minPar.at(0) << "\t" << minPar.at(1) << std::endl;
    double intens = pPIF_->intensity(masssq, minPar.at(0), minPar.at(1));

    chi += fabs(masssq - intens)/(double)nEvents; //TODO: real Chi2
  }

  return chi;
}
