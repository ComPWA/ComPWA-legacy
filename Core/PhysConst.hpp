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

enum class QuantumNumbers {SPIN, ORBITAL_ANGULAR_MOMENTUM, ISOSPIN, CHARGE, PARITY, CPARITY};

struct ParticleProperties {
  std::string name_;
  double mass_;
  double width_;
  int id_;
  int charge_;
  Spin isospin_;
  Spin spin_;
  int parity_;
  int cparity_;

  ParticleProperties() :
      name_(""), mass_(0.0), width_(0.0), id_(0), charge_(0), parity_(0), cparity_(
          0) {
  }

  bool operator()(const ParticleProperties &lhs) const {
    if (charge_ != lhs.charge_)
      return false;
    if (isospin_ != lhs.isospin_)
      return false;
    if (spin_.J_numerator_ != lhs.spin_.J_numerator_)
      return false;
    if (spin_.J_denominator_ != lhs.spin_.J_denominator_)
      return false;
    if (parity_ != lhs.parity_)
      return false;
    if (cparity_ != lhs.cparity_)
      return false;

    return true;
  }
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

  void initQuantumNumberMapping();

  const Constant& findConstant(const std::string& name) const;
  const ParticleProperties& findParticle(const std::string& name) const;
  const ParticleProperties& findParticle(int pid) const;
  std::vector<ParticleProperties> findParticlesWithQN(
      const ParticleProperties& qn) const;
  bool particleExists(const std::string& name) const;

  std::string getQuantumNumberName(
      const QuantumNumbers& qn_type) const;
  QuantumNumbers getQuantumNumberEnum(
      const std::string& qn_name) const;

private:
  PhysConst();

  void readFile();

  std::map<QuantumNumbers, std::string> quantum_number_key_name_mapping_;
  std::map<std::string, QuantumNumbers> name_quantum_number_key_mapping_;

  std::string particleFileName;
  std::string particleDefaultFileName;
  std::string constantFileName;
  std::string constantDefaultFileName;
  std::vector<ParticleProperties> particle_properties_list_;
  std::vector<Constant> constants_list_;
};

} /* namespace ComPWA */

#endif /* PHYSCONST_HPP_ */
