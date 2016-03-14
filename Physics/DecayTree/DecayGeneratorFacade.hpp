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

#ifndef PHYSICS_DECAYTREE_DECAYGENERATORFACADE_HPP_
#define PHYSICS_DECAYTREE_DECAYGENERATORFACADE_HPP_

#include <vector>
#include <string>

#include "Physics/DecayTree/DecayGenerator.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayTree {

class DecayGeneratorFacade {
  DecayGenerator &decay_generator_;

public:
  DecayGeneratorFacade(DecayGenerator &decay_generator);
  virtual ~DecayGeneratorFacade();

  void setAllowedSpinQuantumNumbers(const ComPWA::QuantumNumbers& qn_type,
      const std::vector<unsigned int>& spin_numerators,
      unsigned int spin_denominator,
      const std::vector<ComPWA::QuantumNumbers>& required_qns,
      const QuantumNumberTypes& type = QuantumNumberTypes::SINGLE_PARTICLE_BASED) const;

  void setAllowedIntQuantumNumbers(const ComPWA::QuantumNumbers& qn_type,
      const std::vector<int>& int_qn_values,
      const std::vector<ComPWA::QuantumNumbers>& required_qns,
      const QuantumNumberTypes& type = QuantumNumberTypes::SINGLE_PARTICLE_BASED) const;

  std::vector<std::string> convertQNTypeListToQNStringList(
      const std::vector<ComPWA::QuantumNumbers>& required_qns) const;

  void setConservedQuantumNumbers(
      const std::vector<ComPWA::QuantumNumbers>& conserved_qn) const;

  /* void setAllowedSpins(const std::vector<unsigned int>& spin_numerators,
   unsigned int spin_denominator = 1);
   void setAllowedIsospins(const std::vector<unsigned int>& isospin_numerators,
   unsigned int isospin_denominator = 1);
   void setAllowedCharges(const std::vector<int>& charges);
   void setAllowedParities(const std::vector<int>& parities);
   void setAllowedCParities(const std::vector<int>& cparities);*/
};

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_DECAYTREE_DECAYGENERATORFACADE_HPP_ */
