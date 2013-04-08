//! Internal container representing a particle.
/*! \class PWAParticle
 * @file PWAParticle.hpp
 * This class provides a internal container for information of a particle. The
 * class provides the momentum 4-vector and pid of the particle.
*/

#ifndef _PWAPARTICLE_HPP_
#define _PWAPARTICLE_HPP_

#include <vector>

class PWAParticle
{

public:

  PWAParticle(){
    px = 0;
    py = 0;
    pz = 0;
    E = 0;
    pid = 0;
  }

  PWAParticle(double inPx, double inPy, double inPz, double inE, int inpid=0){
    px = inPx;
    py = inPy;
    pz = inPz;
    E = inE;
    pid = inpid;
  }

  PWAParticle(const PWAParticle& inParticle){
    px = inParticle.px;
    py = inParticle.py;
    pz = inParticle.pz;
    E = inParticle.E;
    pid = inParticle.pid;
  }

  virtual ~PWAParticle(){ /* nothing */	}

  virtual const inline double getPx() const {return px;}
  virtual const inline double getPy() const {return py;}
  virtual const inline double getPz() const {return pz;}
  virtual const inline double getE() const {return E;}
  virtual const inline int getPid() const {return pid;}

protected:
  double px, py, pz, E;
  int pid;
  //TODO: other particle info?

};

#endif
