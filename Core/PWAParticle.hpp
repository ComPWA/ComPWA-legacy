//! Internal container representing a particle.
/*! \class PWAParticle
 * @file PWAParticle.hpp
 * This class provides a internal container for information of a particle. The
 * class provides the momentum 4-vector and pid of the particle.
*/

#ifndef PWAPARTICLE_HPP_
#define PWAPARTICLE_HPP_

#include <vector>

class PWAParticle
{

public:

	PWAParticle()
	  {
		 px = 0;
		 py = 0;
		 pz = 0;
		 E = 0;
		 pid = 0;
	  }

	PWAParticle(double inPx, double inPy, double inPz, double inE, int inpid=0)
	  {
		 px = inPx;
		 py = inPy;
		 pz = inPz;
		 E = inE;
		 pid = inpid;
	  }

	PWAParticle(const PWAParticle& inParticle)
	  {
		 px = inParticle.getPx();
		 py = inParticle.getPy();
		 pz = inParticle.getPz();
		 E = inParticle.getE();
		 pid = inParticle.getPid();
	  }

  virtual ~PWAParticle()
	{ /* nothing */	}

  virtual const inline unsigned int getPx() const {return px;}
  virtual const inline unsigned int getPy() const {return px;}
  virtual const inline unsigned int getPz() const {return px;}
  virtual const inline unsigned int getE() const {return E;}
  virtual const inline unsigned int getPid() const {return pid;}

protected:
  double px, py, pz, E;
  int pid;
  //TODO: other particle info?

};

#endif
