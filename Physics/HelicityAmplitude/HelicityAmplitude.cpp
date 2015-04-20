#include "HelicityAmplitude.hpp"

#include <cmath>

namespace HelicityFormalism {

HelicityAmplitude::HelicityAmplitude() {
	// TODO Auto-generated constructor stub

}

HelicityAmplitude::~HelicityAmplitude() {
	// TODO Auto-generated destructor stub
}

bool HelicityAmplitude::validate() const {
}

std::complex<double> HelicityAmplitude::evaluate() const {
  // those are the two variables that we get from the user
  double theta;
  double phi;
  // now transform those into the cms frame of our initial state

	double spin_factor = sqrt((2*initial_state.J + 1)/(4*M_PI));
	std::complex<double> Djm = Wigner_D(phi, theta, -phi, initial_state.J, initial_state.M, final_state.helicity_particle1-final_state.helicity_particle2);

	return spin_factor * Djm;
}

} /* namespace HelicityFormalism */
