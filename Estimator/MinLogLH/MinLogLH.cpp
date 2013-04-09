#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
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
    Event theEvent(pDIF_->getEvent(evt));
    const Particle &a(theEvent.getParticle(0));
    const Particle &b(theEvent.getParticle(1));
   /* TODO: try read exceptions
    if( !(pDIF_->getEvent(evt, theEvent)) ){
      std::cout << "EIFChiOneD::controlParameter: Event not readable!" << std::endl; //TODO Exception
      continue;
    }*/
    /*if( !(theEvent.getParticle(0,a)) ){
      std::cout << "EIFChiOneD::controlParameter: Particle A not readable!" << std::endl; //TODO Exception
      continue;
    }
    if( !(theEvent.getParticle(1,b)) ){
      std::cout << "EIFChiOneD::controlParameter: Particle B not readable!" << std::endl; //TODO Exception
      continue;
    }*/

    double masssq = 0;
    masssq += (pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2));

    std::vector<double> x;
    x.push_back(sqrt(masssq));
    double intens = pPIF_->intensity(x, minPar);

    lh -= (log(intens/norm));
  }

  return lh;
}
