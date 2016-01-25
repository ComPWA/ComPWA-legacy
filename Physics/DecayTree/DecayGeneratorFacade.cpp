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
#include "Core/PhysConst.hpp"

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

void DecayGeneratorFacade::setAllowedSpinQuantumNumbers(
    const ComPWA::QuantumNumbers& qn_type,
    const std::vector<unsigned int>& spin_numerators,
    unsigned int spin_denominator) const {
  AllowedQuantumNumbers<ComPWA::Spin> spin;

  spin.quantum_number_name_ =
      ComPWA::PhysConst::Instance().getQuantumNumberName(qn_type);
  Spin s;
  s.J_denominator_ = spin_denominator;
  for (unsigned int i = 0; i < spin_numerators.size(); ++i) {
    s.J_numerator_ = spin_numerators[i];
    for (int spinz_num = -spin_numerators[i];
        spinz_num <= (int) spin_numerators[i]; ++spinz_num) {
      s.J_z_numerator_ = spinz_num;
      spin.allowed_values_.push_back(s);
    }
  }
  decay_generator_.allowed_spin_like_quantum_numbers_.push_back(spin);
}

void DecayGeneratorFacade::setAllowedIntQuantumNumbers(
    const ComPWA::QuantumNumbers& qn_type,
    const std::vector<int>& int_qn_values) const {
  AllowedQuantumNumbers<int> int_qn;
  int_qn.quantum_number_name_ =
      ComPWA::PhysConst::Instance().getQuantumNumberName(qn_type);
  int_qn.allowed_values_ = int_qn_values;
  decay_generator_.allowed_integer_like_quantum_numbers_.push_back(int_qn);
}

} /* namespace DecayTree */
} /* namespace ComPWA */
