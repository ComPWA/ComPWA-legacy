//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//    Stefan Pflueger - initial implementation
//-------------------------------------------------------------------------------

#include "Core/Utility.hpp"

#include "Physics/DecayTree/DecayGeneratorFacade.hpp"
#include "Physics/DecayTree/DecayGenerator.hpp"

namespace ComPWA {
namespace DecayTree {

DecayGeneratorFacade::DecayGeneratorFacade(DecayGenerator &decay_generator) :
    decay_generator_(decay_generator) {
}

DecayGeneratorFacade::~DecayGeneratorFacade() {
  // TODO Auto-generated destructor stub
}

void DecayGeneratorFacade::setAllowedSpins(
    const std::vector<unsigned int>& spin_numerators,
    unsigned int spin_denominator) {
  Spin s;
  s.J_denominator_ = spin_denominator;
  for (unsigned int i = 0; i < spin_numerators.size(); ++i) {
    s.J_numerator_ = spin_numerators[i];
    for (int spinz_num = -spin_numerators[i]; spinz_num <= (int)spin_numerators[i];
        ++spinz_num) {
      s.J_z_numerator_ = spinz_num;
      decay_generator_.allowed_spins_.push_back(s);
    }
  }
}

void DecayGeneratorFacade::setAllowedIsospins(
    const std::vector<unsigned int>& isospin_numerators,
    unsigned int isospin_denominator) {
  Spin s;
  s.J_denominator_ = isospin_denominator;
  for (unsigned int i = 0; i < isospin_numerators.size(); ++i) {
    s.J_numerator_ = isospin_numerators[i];
    for (int spinz_num = -isospin_numerators[i];
        spinz_num <= (int)isospin_numerators[i]; ++spinz_num) {
      s.J_z_numerator_ = spinz_num;
      decay_generator_.allowed_isospins_.push_back(s);
    }
  }
}

void DecayGeneratorFacade::setAllowedCharges(
    const std::vector<int>& charges) {
  decay_generator_.allowed_charges_ = charges;
}
void DecayGeneratorFacade::setAllowedParities(
    const std::vector<int>& parities) {
  decay_generator_.allowed_parities_ = parities;
}
void DecayGeneratorFacade::setAllowedCParities(
    const std::vector<int>& cparities) {
  decay_generator_.allowed_cparites_ = cparities;
}

} /* namespace DecayTree */
} /* namespace ComPWA */
