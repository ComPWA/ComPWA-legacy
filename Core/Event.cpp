#include <vector>
#include <string>
#include "Core/Particle.hpp"
#include "Core/Event.hpp"

Event::Event():fWeight(1.),fName(""){

}

Event::Event(const std::string& name):fWeight(1.),fName(name){

}

Event::Event(const double inWeight, const std::string& name=""):fWeight(inWeight),fName(name){

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
