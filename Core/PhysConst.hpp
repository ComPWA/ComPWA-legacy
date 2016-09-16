//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

//! PhysConst provides particle information to the framework
/*!
 * @file PhysConst.hpp
 *\class PhysConst
 *      PhysConst reads particle properties from an xml file. And provides
 *      it to the framework via a singleton structure. The structure of
 *      the xml file is defined in particleSchema.xsd. It is also intended
 *      to store other physical constants in this way.
 */

#ifndef PHYSCONST_HPP_
#define PHYSCONST_HPP_

#include <iostream>
#include <vector>

#include <Core/Parameter.hpp>
#include <Core/Utility.hpp>

namespace ComPWA {

enum class QuantumNumberIDs {
  SPIN, ORBITAL_ANGULAR_MOMENTUM, ISOSPIN, CHARGE, PARITY, CPARITY
};
enum class QuantumNumberType {
  SPIN_LIKE, INTEGER_LIKE, DOUBLE_LIKE
};

class QuantumNumberTranslator {
  std::map<QuantumNumberIDs, std::string> quantum_number_key_name_mapping_;
  std::map<std::string, QuantumNumberIDs> name_quantum_number_key_mapping_;
  std::map<QuantumNumberIDs, QuantumNumberType> quantum_number_key_type_mapping_;

  QuantumNumberTranslator();
public:
  static QuantumNumberTranslator& Instance() {
    static QuantumNumberTranslator instance;
    return instance;
  }
  ~QuantumNumberTranslator();

  QuantumNumberTranslator(QuantumNumberTranslator const&) = delete;
  void operator=(QuantumNumberTranslator const&) = delete;

  QuantumNumberType getQuantumNumberType(const std::string& qn_name) const;
  std::string getQuantumNumberName(const QuantumNumberIDs& qn_type) const;
  QuantumNumberIDs getQuantumNumberEnum(const std::string& qn_name) const;
};

struct SpinWave {
  std::map<std::string, ComPWA::Spin> spin_like_quantum_numbers_;
  std::map<std::string, int> integer_like_quantum_numbers_;
  std::map<std::string, double> double_like_quantum_numbers_;

  ComPWA::Spin getSpinLikeQuantumNumber(QuantumNumberIDs qn_id) const;
  int getIntLikeQuantumNumber(QuantumNumberIDs qn_id) const;
  double getDoubleLikeQuantumNumber(QuantumNumberIDs qn_id) const;

  bool operator==(const SpinWave& rhs) const {
    if (spin_like_quantum_numbers_ != rhs.spin_like_quantum_numbers_)
      return false;
    if (integer_like_quantum_numbers_ != rhs.integer_like_quantum_numbers_)
      return false;
    if (double_like_quantum_numbers_ != rhs.double_like_quantum_numbers_)
      return false;
    return true;
  }
  // lets define a modulus operator for our special comparison
  // we want to check if all the properties of rhs are available in this
  bool supersetOf(const SpinWave& rhs) const {
    for (auto const& spin_like_qn : rhs.spin_like_quantum_numbers_) {
      bool nothing_found(true);
      auto result = spin_like_quantum_numbers_.find(spin_like_qn.first);
      if (result != spin_like_quantum_numbers_.end()) {
        if (1.0 * result->second.GetNumerator() / result->second.GetDenominator()
            == 1.0 * spin_like_qn.second.GetNumerator()
                / spin_like_qn.second.GetDenominator()) {
          if (result->second.UseZ()) {
            if (1.0 * result->second.GetZNumerator()
                / result->second.GetDenominator()
                == 1.0 * spin_like_qn.second.GetZNumerator()
                    / spin_like_qn.second.GetDenominator()) {
              nothing_found = false;
            }
          }
          else
            nothing_found = false;
        }
      }
      if (nothing_found)
        return false;
    }
    for (auto const& int_like_qn : rhs.integer_like_quantum_numbers_) {
      auto result = integer_like_quantum_numbers_.find(int_like_qn.first);
      if (result != integer_like_quantum_numbers_.end()) {
        if (result->second != int_like_qn.second) {
          return false;
        }
      }
    }
    for (auto const& double_like_qn : rhs.double_like_quantum_numbers_) {
      auto result = double_like_quantum_numbers_.find(double_like_qn.first);
      if (result != double_like_quantum_numbers_.end()) {
        if (result->second != double_like_qn.second) {
          return false;
        }
      }
    }

    return true;
  }

  friend std::ostream& operator<<(std::ostream& stream, const SpinWave& sw) {
    for (auto iter = sw.spin_like_quantum_numbers_.begin();
        iter != sw.spin_like_quantum_numbers_.end(); ++iter) {
      stream << iter->first << ": " << iter->second.GetNumerator()<< "/"
          << iter->second.GetDenominator()<< " (z="
          << iter->second.GetZNumerator()<< ")\n";
    }
    for (auto iter = sw.integer_like_quantum_numbers_.begin();
        iter != sw.integer_like_quantum_numbers_.end(); ++iter) {
      stream << iter->first << ": " << iter->second << "\n";
    }
    return stream;
  }
};

struct ParticleProperties: public SpinWave {
  std::string name_;
  double mass_;
  double width_;
  int id_;

  ParticleProperties() :
      name_(""), mass_(0.0), width_(0.0), id_(0) {
  }

  /*bool operator()(const ParticleProperties &lhs) const {
   if (spin_like_quantum_numbers_ != lhs.spin_like_quantum_numbers_)
   return false;
   if (integer_like_quantum_numbers_ != lhs.integer_like_quantum_numbers_)
   return false;
   if (double_like_quantum_numbers_ != lhs.double_like_quantum_numbers_)
   return false;

   return true;
   }*/
};

struct Constant {
  std::string name_;
  double value_;
  double error_;

  Constant() :
      name_(""), value_(0.0), error_(0.0) {
  }
};

class PhysConst {
public:
  //! returns existing instance of PhysConst, or creates a new one
  static PhysConst& Instance() {
    static PhysConst instance;
    return instance;
  }
  ~PhysConst();

  PhysConst(PhysConst const&) = delete;
  void operator=(PhysConst const&) = delete;

  const Constant& findConstant(const std::string& name) const;
  const ParticleProperties& findParticle(const std::string& name) const;
  const ParticleProperties& findParticle(int pid) const;
  std::vector<ParticleProperties> findParticlesWithQN(
      const ParticleProperties& qn) const;
  bool particleExists(const std::string& name) const;

private:
  PhysConst();

  void readFile();

  std::string particleFileName;
  std::string particleDefaultFileName;
  std::string constantFileName;
  std::string constantDefaultFileName;
  std::vector<ParticleProperties> particle_properties_list_;
  std::vector<Constant> constants_list_;
};

} /* namespace ComPWA */

#endif /* PHYSCONST_HPP_ */
