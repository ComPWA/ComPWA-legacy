//! Internal container for event information.
/*! \class Event
 * @file Event.hpp
 * This class provides a internal container for event-based information. The
 * class provides a list of particles of the event.
*/

#ifndef _Event_HPP_
#define _Event_HPP_

#include <vector>
#include "Core/Particle.hpp"

class Event{

public:

  Event();

  Event(const double inWeight);

  virtual void addParticle(Particle inParticle);

  virtual ~Event();

  virtual const inline unsigned int getNParticles() { return fParticles.size(); }
  virtual const Particle& getParticle(const unsigned int id);

protected:
  std::vector<Particle> fParticles;
  double fWeight;
  //Particle fParticleB;
  //TODO: other event info?

};

#endif
