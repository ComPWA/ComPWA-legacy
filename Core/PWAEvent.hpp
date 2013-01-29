//! Internal container for event information.
/*! \class PWAEvent
 * @file PWAEvent.hpp
 * This class provides a internal container for event-based information. The
 * class provides a list of particles of the event.
*/

#ifndef _PWAEvent_HPP_
#define _PWAEvent_HPP_

#include <vector>
#include "Core/PWAParticle.hpp"

class PWAEvent{

public:

  PWAEvent();

  PWAEvent(const double inWeight);

  virtual void addParticle(PWAParticle inParticle);

  virtual ~PWAEvent();

  virtual const inline unsigned int getNParticles() { return fParticles.size(); }
  virtual const int getParticle(const unsigned int id, PWAParticle& out);

protected:
  std::vector<PWAParticle> fParticles;
  double fWeight;
  //PWAParticle fParticleB;
  //TODO: other event info?

};

#endif
