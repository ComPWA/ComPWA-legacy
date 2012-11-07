//! Internal container for event information.
/*! \class PWAEvent
 * @file PWAEvent.hpp
 * This class provides a internal container for event-based information. The
 * class provides a list of particles of the event.
*/

#ifndef PWAEvent_HPP_
#define PWAEvent_HPP_

#include <vector>
#include <PWAParticle.hpp>

class PWAEvent{

public:

  PWAEvent(){
  }

  virtual void addParticle(const PWAParticle inParticle){
    fParticles.push_back(inParticle);
  }

  virtual ~PWAEvent() { /* nothing */	}

  virtual const inline unsigned int getNParticles() { return fParticles.size(); }
  virtual const int getParticle(const unsigned int id, PWAParticle& out){
    if(id>=getNParticles()){
      //TODO Exception
      return 0;
    } else {
      out = fParticles.at(id);
      return 1;
    }
    return 0;
  }

protected:
  std::vector<PWAParticle> fParticles;
  //PWAParticle fParticleB;
  //TODO: other event info?

};

#endif
