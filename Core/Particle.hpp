//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Internal container representing a particle.
/*! \class Particle
 * @file Particle.hpp
 * This class provides a internal container for information of a particle. The
 * class provides the momentum 4-vector and pid of the particle.
 */

#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace ComPWA {
class Particle
{
public:
	//! Default constructor
	Particle(double inPx=0, double inPy=0, double inPz=0,
			double inE=0, int inpid=0, int c=0);

	//! Default destructor
	virtual ~Particle();

	//! Invariant mass of this particle and @param inParticle
	inline double invariantMass(const Particle& in) const
	{
		return ( std::pow(getE()+in.getE(),2)
		-std::pow(getPx()+in.getPx() ,2)
		-std::pow(getPy()+in.getPy() ,2)
		-std::pow(getPz()+in.getPz() ,2)
		);
	}
	/**! Invariant mass
	 *
	 * @param inPa Particle A
	 * @param inPb Particle B
	 * @return
	 */
	static double invariantMass(const Particle& inPa, const Particle& inPb);

	//!Get magnitude of three momentum
	double getThreeMomentum() const { return sqrt(px*px+py*py+pz*pz); }
	//!Get invariant mass
	double getInvMass() const { return sqrt(E*E-px*px+py*py+pz*pz); }

	virtual const inline void setPx(double _px) {px=_px;}
	virtual const inline void setPy(double _py) {py=_py;}
	virtual const inline void setPz(double _pz) {pz=_pz;}
	virtual const inline void setE(double _e) {E=_e;}
	virtual const inline void setPid(int _pid) {pid=_pid;}
	virtual const inline void setCharge(int _c) { charge=_c; }
	virtual const inline double getPx() const {return px;}
	virtual const inline double getPy() const {return py;}
	virtual const inline double getPz() const {return pz;}
	virtual const inline double getE() const {return E;}
	virtual const inline int getPid() const {return pid;}
	virtual const inline int getCharge() const {return charge;}

	friend std::ostream& operator<< (std::ostream& stream, const Particle& p);

	Particle& operator+=(const Particle& rhs) ;

	Particle operator+(const Particle& c2) ;

	inline double getMassSquare() const {
		return E*E-px*px-py*py-pz*pz;
	}

	//protected:
	double px, py, pz, E;
	int pid;
	int charge;
};

} /* namespace ComPWA */
#endif
