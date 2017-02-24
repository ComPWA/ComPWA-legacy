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

#include <boost/property_tree/ptree.hpp>

#include <Core/Parameter.hpp>
#include <Core/Spin.hpp>

namespace ComPWA {

enum class QuantumNumberIDs {
  SPIN,
  ORBITAL_ANGULAR_MOMENTUM,
  ISOSPIN,
  CHARGE,
  PARITY,
  CPARITY
};
enum class QuantumNumberType { SPIN_LIKE, INTEGER_LIKE, DOUBLE_LIKE };

class QuantumNumberTranslator {
  std::map<QuantumNumberIDs, std::string> quantum_number_key_name_mapping_;
  std::map<std::string, QuantumNumberIDs> name_quantum_number_key_mapping_;
  std::map<QuantumNumberIDs, QuantumNumberType>
      quantum_number_key_type_mapping_;

  QuantumNumberTranslator();

public:
  static QuantumNumberTranslator &Instance() {
    static QuantumNumberTranslator instance;
    return instance;
  }
  ~QuantumNumberTranslator();

  QuantumNumberTranslator(QuantumNumberTranslator const &) = delete;
  void operator=(QuantumNumberTranslator const &) = delete;

  QuantumNumberType getQuantumNumberType(const std::string &qn_name) const;
  std::string getQuantumNumberName(const QuantumNumberIDs &qn_type) const;
  QuantumNumberIDs getQuantumNumberEnum(const std::string &qn_name) const;
};

class SpinWave {
public:
  ComPWA::Spin getSpinLikeQuantumNumber(QuantumNumberIDs qn_id) const;
  int getIntLikeQuantumNumber(QuantumNumberIDs qn_id) const;
  double getDoubleLikeQuantumNumber(QuantumNumberIDs qn_id) const;

  bool operator==(const SpinWave &rhs) const {
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
  bool supersetOf(const SpinWave &rhs) const {
    for (auto const &spin_like_qn : rhs.spin_like_quantum_numbers_) {
      bool nothing_found(true);
      auto result = spin_like_quantum_numbers_.find(spin_like_qn.first);
      if (result != spin_like_quantum_numbers_.end()) {
        if (1.0 * result->second.GetNumerator() /
                result->second.GetDenominator() ==
            1.0 * spin_like_qn.second.GetNumerator() /
                spin_like_qn.second.GetDenominator()) {
          if (result->second.UseZ()) {
            if (1.0 * result->second.GetZNumerator() /
                    result->second.GetDenominator() ==
                1.0 * spin_like_qn.second.GetZNumerator() /
                    spin_like_qn.second.GetDenominator()) {
              nothing_found = false;
            }
          } else
            nothing_found = false;
        }
      }
      if (nothing_found)
        return false;
    }
    for (auto const &int_like_qn : rhs.integer_like_quantum_numbers_) {
      auto result = integer_like_quantum_numbers_.find(int_like_qn.first);
      if (result != integer_like_quantum_numbers_.end()) {
        if (result->second != int_like_qn.second) {
          return false;
        }
      }
    }
    for (auto const &double_like_qn : rhs.double_like_quantum_numbers_) {
      auto result = double_like_quantum_numbers_.find(double_like_qn.first);
      if (result != double_like_quantum_numbers_.end()) {
        if (result->second != double_like_qn.second) {
          return false;
        }
      }
    }

    return true;
  }

  friend std::ostream &operator<<(std::ostream &stream, const SpinWave &sw) {
    for (auto iter = sw.spin_like_quantum_numbers_.begin();
         iter != sw.spin_like_quantum_numbers_.end(); ++iter) {
      stream << iter->first << ": " << iter->second.GetNumerator() << "/"
             << iter->second.GetDenominator()
             << " (z=" << iter->second.GetZNumerator() << ")\n";
    }
    for (auto iter = sw.integer_like_quantum_numbers_.begin();
         iter != sw.integer_like_quantum_numbers_.end(); ++iter) {
      stream << iter->first << ": " << iter->second << "\n";
    }
    return stream;
  }

protected:
  std::map<std::string, ComPWA::Spin> spin_like_quantum_numbers_;
  std::map<std::string, int> integer_like_quantum_numbers_;
  std::map<std::string, double> double_like_quantum_numbers_;
};

class Properties {
public:
  Properties(std::string name = "test", int id = -999) : _name(name), _id(id){};

  void SetName(std::string n) { _name = n; }
  std::string GetName() const { return _name; }

  void SetId(int id) { _id = id; }
  int GetId() const { return _id; }

protected:
  std::string _name;
  int _id;
};

class ParticleProperties : public Properties {
public:
  ParticleProperties(std::string name = "test", int id = -999)
      : Properties(name, id){};
  
  ParticleProperties(boost::property_tree::ptree pt){
    SetName( pt.get<std::string>("<xmlattr>.Name") );
    SetId( pt.get<int>("Id") );
    SetMass( DoubleParameterFactory( pt.get_child("Mass") ) );
    SetSpin( ComPWA::Spin(pt.get<double>("Spin")) );
    SetIsoSpin( ComPWA::Spin(pt.get<double>("IsoSpin")) );
    SetIsoSpinZ( ComPWA::Spin(pt.get<double>("IsoSpinZ")) );
    SetParity( pt.get<double>("Parity") );
    SetCparity( pt.get<double>("Cparity") );
    SetGparity( pt.get<double>("Gparity") );
  }

  void SetMass(double m) { _mass.SetValue(m); }
  double GetMass() const { return _mass.GetValue(); }

  void SetMassPar(ComPWA::DoubleParameter m) { _mass = m; }
  ComPWA::DoubleParameter GetMassPar() const { return _mass; }

  void SetSpin(ComPWA::Spin s) { _spin = s; }
  ComPWA::Spin GetSpin() const { return _spin; }

  void SetIsoSpin(ComPWA::Spin s) { _isoSpin = s; }
  ComPWA::Spin GetIsoSpin() const { return _isoSpin; }

  void SetIsoSpinZ(ComPWA::Spin s) { _isoSpinZ = s; }
  ComPWA::Spin GetIsoSpinZ() const { return _isoSpinZ; }

  void SetParity(int p) { _parity = p; }
  int GetParity() const { return _parity; }

  void SetCparity(int p) { _cparity = p; }
  int GetCparity() const { return _cparity; }

  void SetGparity(int p) { _gparity = p; }
  int GetGparity() const { return _gparity; }

protected:
  ComPWA::DoubleParameter _mass;
  ComPWA::Spin _spin;
  int _parity;
  int _cparity;
  int _gparity;
  ComPWA::Spin _isoSpin;
  ComPWA::Spin _isoSpinZ;

  /* Store decay info in property_tree. The tree is later on passed to the
   * respective class. */
  boost::property_tree::ptree decayInfo;
};

class Constant : public Properties {
public:
  Constant(std::string n="", double value=0.0) : Properties(n), _value(n,value) { }
  Constant(boost::property_tree::ptree pt){};

  void SetValue(double m) { _value.SetValue(m); }
  double GetValue() const { return _value.GetValue(); }

  void SetValuePar(ComPWA::DoubleParameter m) { _value = m; }
  ComPWA::DoubleParameter GetValuePar() const { return _value; }

protected:
  ComPWA::DoubleParameter _value;
};

class PhysConst {
public:
  static PhysConst *CreateInstance(std::string file = "./particles.xml") {
    if (_inst)
        throw std::runtime_error("PhysConst::CreateInstance() | Instance already exists. Use Instance() to access it!");
    _inst = new PhysConst(file);
    return _inst;
  }
  
  static PhysConst *CreateInstance(boost::property_tree::ptree pt) {
    if (_inst)
      throw std::runtime_error("PhysConst::CreateInstance() | Instance already exists. Use Instance() to access it!");
    
    _inst = new PhysConst(pt);
    return _inst;
  }


  static PhysConst *Instance() {
    if (!_inst) {
      throw std::runtime_error("No instance of PhysConst created! "
                               "Create one first!");
    }
    return _inst;
  }

  ~PhysConst();

  
  PhysConst(PhysConst const &) = delete;

  void operator=(PhysConst const &) = delete;

  const Constant &FindConstant(const std::string &name) const;

  const ParticleProperties &FindParticle(const std::string &name) const;

  const ParticleProperties &FindParticle(int pid) const;

//  std::vector<ParticleProperties>
//  findParticlesWithQN(const ParticleProperties &qn) const;

  bool ParticleExists(const std::string &name) const;

protected:
  PhysConst(std::string filePath = "./particles.xml");

  PhysConst(boost::property_tree::ptree pt);
  
  //! Singleton stuff
  static PhysConst *_inst;
  
  void readTree( boost::property_tree::ptree pt );
  
  //! Input file name
  std::string _particleFileName;
  //! Input file name
  std::string _constantFileName;
  
  std::vector<ParticleProperties> _partList;
  std::vector<Constant> _constList;
};

} /* namespace ComPWA */

#endif /* PHYSCONST_HPP_ */
