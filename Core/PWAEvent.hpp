#ifndef PWAEvent_HPP_
#define PWAEvent_HPP_

#include <vector>
#include <PWAParticle.hpp>

using namespace std;

class PWAEvent
{

public:

	PWAEvent()
	  {
	  }

	PWAEvent(PWAParticle inA, PWAParticle inB)
	  {
		fParticleA = inA;
		fParticleB = inB;
	  }

  virtual ~PWAEvent()
	{ /* nothing */	}

  virtual const inline unsigned int getNParticles() {return 2;/*fParticles.size();*/}
  virtual const PWAParticle& getParticle(const unsigned int id)
    {
	  //if(id<0 || id>=getNParticles()){
		  //TODO Exception
		//  return PWAParticle();
	  //}
	  if(id==0) return fParticleA;
	  if(id==1) return fParticleB;
	  return PWAParticle();
    }

protected:
  PWAParticle fParticleA;
  PWAParticle fParticleB;
  //TODO: other event info?

};

#endif
