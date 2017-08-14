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

/// Particle ID.
/// Usually the pid's from PDG are used here:
/// http://pdg.lbl.gov/mc_particle_id_contents.html
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

///
/// \class Constant
/// Represents physical constants.
///
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

///
/// \class PartInfoShort
/// In many cases only very basic information on particles is needed.
/// Therefore, we use this class to store {id, name, mass}
///
class PartInfoShort : public Properties {
public:
  PartInfoShort(std::string name = "test", pid id = -999, double mass = 0)
      : Properties(name, id), _mass(mass){};

  void SetMass(double m) { _mass = m; }
  double GetMass() const { return _mass; }

protected:
  double _mass;
};

///
/// \class ParticleProperties
///
class ParticleProperties : public Properties {
public:
  ParticleProperties(std::string name = "test", pid id = -999)
      : Properties(name, id){};

  ParticleProperties(boost::property_tree::ptree pt);
  
  virtual boost::property_tree::ptree Save();

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

  /// Store decay info in property_tree. The tree is later on passed to the
  /// respective class.
  boost::property_tree::ptree _decayInfo;
};

/// A map of particle properties is used everywhere where particle information
/// is needed. Properties are accessed by the particle name.
/// Note: Propably would be better to access particles by their pid?
typedef std::map<std::string, ParticleProperties> PartList;

/// Search particle \p list for a specific particle \p id.
/// The first entry in the list is returned. Be careful in case that multiple
/// particles have the same pid.
// const ParticleProperties &FindParticle(std::shared_ptr<PartList> list, pid
// id);
inline const ParticleProperties &FindParticle(std::shared_ptr<PartList> list,
                                              pid id) {

  // position in map
  int result = -1;
  for (int i = 0; i < list->size(); ++i) {
    // if a match is found we skip the rest
    if (result >= 0)
      continue;

    auto it = list->begin();
    std::advance(it, i);
    if (it->second.GetId() == id) {
      result = i;
    }
  }
  if (result < 0) {
    std::stringstream ss;
    ss << "FindParticle() | Particle id=" << id << " not found in list!";
    throw std::runtime_error(ss.str());
  }
  auto r = list->begin();
  std::advance(r, result);

  return r->second;
}

/// Read list of particles from a boost::property_tree
inline void ReadParticles(std::shared_ptr<PartList> list,
                          boost::property_tree::ptree pt) {

  auto particleTree = pt.get_child_optional("ParticleList");
  if (!particleTree)
    return;

  for (auto const &v : particleTree.get()) {
    auto tmp = ParticleProperties(v.second);
    auto p = std::make_pair(tmp.GetName(), tmp);
    auto last = list->insert(p);

    if (!last.second) {
      LOG(info) << "ReadParticles() | Particle " << last.first->first
                << "already exists in list. We overwrite its parameters!";
      last.first->second = tmp;
    }
    tmp = last.first->second;

    // cparity is optional
    double cparity = 0.0;
    try {
      cparity = tmp.GetQuantumNumber("Cparity");
    } catch (std::exception &ex) {
    }

    LOG(info) << "ReadParticles() | Particle " << tmp.GetName()
              << " (id=" << tmp.GetId() << ") "
              << " J(PC)=" << tmp.GetSpinQuantumNumber("Spin") << "("
              << tmp.GetQuantumNumber("Parity") << cparity << ") "
              << " mass=" << tmp.GetMass()
              << " decayType=" << tmp.GetDecayType();
  }

  return;
}

/// Save particle list to boost::property_tree
inline boost::property_tree::ptree
SaveParticles(std::shared_ptr<PartList> list) {
  boost::property_tree::ptree pt;
  
  for( auto& i: *list.get() ){
    pt.add_child("Particle",i.second.Save());
  }
  boost::property_tree::ptree pt2;
  pt2.add_child("ParticleList",pt);
  return pt2;
}

} // Namespace ComPWA

#endif
