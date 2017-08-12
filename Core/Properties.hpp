// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef Properties_h
#define Properties_h

#include <vector>

#include <boost/property_tree/ptree.hpp>

#include <Core/Exceptions.hpp>
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

  ParticleProperties(boost::property_tree::ptree pt);

  double GetMass() const { return _mass.GetValue(); }

  ComPWA::DoubleParameter GetMassPar() const { return _mass; }

  int GetQuantumNumber(std::string type) const;

  ComPWA::Spin GetSpinQuantumNumber(std::string type) const;

  boost::property_tree::ptree GetDecayInfo() const { return _decayInfo; }

  std::string GetDecayType() const {
    return _decayInfo.get<std::string>("<xmlattr>.Type");
  }

protected:
  ComPWA::DoubleParameter _mass;
  std::map<std::string, int> intQuantumNumbers_;
  std::map<std::string, ComPWA::Spin> spinQuantumNumbers_;

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

} // Namespace ComPWA
#endif /* Properties_h */
