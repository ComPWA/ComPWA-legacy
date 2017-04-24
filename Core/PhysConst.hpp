
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

  /*! Particle ID.
   * Usually the pid's from PDG are used here: 
   * http://pdg.lbl.gov/mc_particle_id_contents.html
   */
  typedef int pid;
  
class Properties {
public:
  Properties(std::string name = "test", pid id = -999) : _name(name), _id(id){};

  void SetName(std::string n) { _name = n; }
  std::string GetName() const { return _name; }

  void SetId(pid id) { _id = id; }
  pid GetId() const { return _id; }

protected:
  std::string _name;
  pid _id;
};

class ParticleProperties : public Properties {
public:
  ParticleProperties(std::string name = "test", pid id = -999)
      : Properties(name, id){};

  ParticleProperties(boost::property_tree::ptree pt)
      : Properties(pt.get<std::string>("<xmlattr>.Name"), pt.get<pid>("Id")) {

    SetMassPar(DoubleParameterFactory(pt.get_child("Mass")));
    // Parameters optional. Shall we require them?
    SetSpin(ComPWA::Spin(pt.get<double>("Spin", 0.)));
    SetIsoSpin(ComPWA::Spin(pt.get<double>("IsoSpin", 0.)));
    SetIsoSpinZ(ComPWA::Spin(pt.get<double>("IsoSpinZ", 0.)));
    SetParity(pt.get<double>("Parity", 0));
    SetCparity(pt.get<double>("Cparity", 0));
    SetGparity(pt.get<double>("Gparity", 0));

    auto decayInfo = pt.get_child_optional("DecayInfo");
    if (decayInfo)
      SetDecayInfo(decayInfo.get());
    else {
      _decayInfo.put("<xmlattr>.Type", "stable");
    }
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

  void SetDecayInfo(boost::property_tree::ptree pt) { _decayInfo = pt; }
  boost::property_tree::ptree GetDecayInfo() const { return _decayInfo; }

  std::string GetDecayType() const {
    return _decayInfo.get<std::string>("<xmlattr>.Type");
  }

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
  boost::property_tree::ptree _decayInfo;
};

class Constant : public Properties {
public:
  Constant(std::string n = "", double value = 0.0)
      : Properties(n), _value(n, value) {}
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
      throw std::runtime_error("PhysConst::CreateInstance() | Instance already "
                               "exists. Use Instance() to access it!");
    _inst = new PhysConst(file);
    return _inst;
  }

  static PhysConst *CreateInstance(boost::property_tree::ptree pt) {
    if (_inst)
      throw std::runtime_error("PhysConst::CreateInstance() | Instance already "
                               "exists. Use Instance() to access it!");

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

  const ParticleProperties &FindParticle(pid id) const;

  bool ParticleExists(const std::string &name) const;

protected:
  PhysConst(std::string filePath = "./particles.xml");

  PhysConst(boost::property_tree::ptree pt);

  //! Singleton stuff
  static PhysConst *_inst;

  void readTree(boost::property_tree::ptree pt);

  //! Input file name
  std::string _particleFileName;
  //! Input file name
  std::string _constantFileName;

  std::vector<ParticleProperties> _partList;
  std::vector<Constant> _constList;
};

} /* namespace ComPWA */

#endif /* PHYSCONST_HPP_ */
