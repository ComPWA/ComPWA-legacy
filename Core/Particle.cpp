/*
 * Particle.cpp
 *
 *  Created on: Nov 16, 2015
 *      Author: weidenka
 */


#include "Core/Particle.hpp"


std::ostream& operator<< (std::ostream& stream, const Particle& p){
	stream<< "Particle id="<<p.pid<<" charge="<<p.charge
			<<" p4=("<<p.px<<","<<p.py<<","<<p.pz<<","<<p.px<<")";
	return stream;
}
