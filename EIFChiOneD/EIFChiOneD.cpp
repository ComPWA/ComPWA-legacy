#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "EIFChiOneD.hpp"
#include "PWAEvent.hpp"
#include "PWAParticle.hpp"
#include "PWAParameter.hpp"

EIFChiOneD::EIFChiOneD(std::shared_ptr<PIFBase> inPIF, std::shared_ptr<DIFBase> inDIF) : pPIF_(inPIF), pDIF_(inDIF){

}

EIFChiOneD::~EIFChiOneD(){

}

double EIFChiOneD::controlParameter(std::vector<PWAParameter<double> >& minPar){
  unsigned int nEvents = pDIF_->getNEvents();

 // std::cout << std::endl << "Test" << std::endl;
 // std::cout << "Events: \t" << nEvents << std::endl;
 // std::cout << "Parameter" << std::endl;
 // for(unsigned int i=0; i<minPar.size(); i++)
 //   std::cout << i << "-t Par " << minPar.at(i) << std::endl;

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

    double masssq = 0;
    masssq += pow(a.getE()+b.getE(),2);
    masssq -= pow(a.getPx()+b.getPx() ,2);
    masssq -= pow(a.getPy()+b.getPy() ,2);
    masssq -= pow(a.getPz()+b.getPz() ,2);

    //std::cout << "Event, Par0, Par1: \t" << evt << "\t" << minPar.at(0) << "\t" << minPar.at(1) << std::endl;
    std::vector<double> x;
    x.push_back(sqrt(masssq));
    double intens = pPIF_->intensity(x, minPar);

    chi += fabs(masssq - intens)/(double)nEvents; //TODO: real Chi2
  }

  return chi;
}
