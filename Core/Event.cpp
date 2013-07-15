#include <vector>
#include "Core/Particle.hpp"
#include "Core/Event.hpp"


Event::Event():fWeight(1.){

}

Event::Event(const double inWeight):fWeight(inWeight){

}

void Event::addParticle(Particle inParticle){
  fParticles.push_back(inParticle);
}

Event::~Event() { /* nothing */	}

const Particle& Event::getParticle(const unsigned int id){
  if(id>=getNParticles()){
    //TODO Exception
    return Particle();
  }
  return fParticles.at(id);
}
