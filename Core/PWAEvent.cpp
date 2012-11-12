#include <vector>
#include "PWAParticle.hpp"
#include "PWAEvent.hpp"


PWAEvent::PWAEvent(){
}

void PWAEvent::addParticle(PWAParticle inParticle){
  fParticles.push_back(inParticle);
}

PWAEvent::~PWAEvent() { /* nothing */	}

const int PWAEvent::getParticle(const unsigned int id, PWAParticle& out){
  if(id>=getNParticles()){
    //TODO Exception
    return 0;
  } else {
    out = fParticles.at(id);
    return 1;
  }
  return 0;
}
