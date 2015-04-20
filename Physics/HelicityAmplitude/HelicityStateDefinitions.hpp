#ifndef HELICITYSTATEDEFINITIONS_HPP_
#define HELICITYSTATEDEFINITIONS_HPP_

#include <qft++.h>

namespace HelicityFormalism {

struct SphericalWaveTwoParticleState {
	LS ls;
	Spin J;
	Spin M;
};

struct PlaneWaveTwoParticleState {
	Spin spin_particle1;
	Spin spin_particle2;
	Spin helicity_particle1;
	Spin helicity_particle2;
};

}

#endif /* HELICITYSTATEDEFINITIONS_HPP_ */
