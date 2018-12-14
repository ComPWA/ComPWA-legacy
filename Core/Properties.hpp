// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef Properties_h
#define Properties_h

#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <Core/Exceptions.hpp>
#include <Core/FitParameter.hpp>
#include <Core/ParameterList.hpp>
#include <Core/Spin.hpp>

namespace ComPWA {

/// Particle ID.
/// Usually the pid's from PDG are used here:
/// http://pdg.lbl.gov/mc_particleId_contents.html
typedef int pid;

class Properties {
public:
  Properties(std::string name = "test", pid id = -999) : Name(name), Id(id){};

  void setName(std::string n) { Name = n; }
  std::string name() const { return Name; }

  void SetId(pid id) { Id = id; }
  pid GetId() const { return Id; }

protected:
  std::string Name;
  pid Id;
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

  void SetValue(double m) { _value.setValue(m); }

  double value() const { return _value.value(); }

  void SetValuePar(ComPWA::FitParameter m) { _value = m; }

  ComPWA::FitParameter GetValuePar() const { return _value; }

protected:
  ComPWA::FitParameter _value;
};

///
/// \class PartInfoShort
/// In many cases only very basic information on particles is needed.
/// Therefore, we use this class to store {id, name, mass}
///
class PartInfoShort : public Properties {
public:
  PartInfoShort(std::string name = "test", pid id = -999, double mass = 0)
      : Properties(name, id), Mass(mass){};

  void SetMass(double m) { Mass = m; }
  double GetMass() const { return Mass; }

protected:
  double Mass;
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

  double GetMass() const { return Mass.value(); }

  ComPWA::FitParameter GetMassPar() const { return Mass; }

  int GetQuantumNumber(std::string type) const;

  ComPWA::Spin GetSpinQuantumNumber(std::string type) const;

  boost::property_tree::ptree GetDecayInfo() const { return DecayInfo; }

  std::string GetDecayType() const {
    return DecayInfo.get<std::string>("<xmlattr>.Type");
  }

protected:
  ComPWA::FitParameter Mass;
  std::map<std::string, int> intQuantumNumbers_;
  std::map<std::string, ComPWA::Spin> spinQuantumNumbers_;

  /// Store decay info in property_tree. The tree is later on passed to the
  /// respective class.
  boost::property_tree::ptree DecayInfo;
};

/// A map of particle properties is used everywhere where particle information
/// is needed. Properties are accessed by the particle name.
/// Note: Propably would be better to access particles by their pid?
typedef std::map<std::string, ParticleProperties> PartList;

inline std::ostream &operator<<(std::ostream &os, const PartList &p) {
  for (auto i : p)
    os << i.first << " [ " << i.second.GetId()
       << " ]: mass = " << i.second.GetMass() << std::endl;
  return os;
}

/// Search particle \p list for a specific particle \p id.
/// The first entry in the list is returned. Be careful in case that multiple
/// particles have the same pid.
// const ParticleProperties &FindParticle(std::shared_ptr<PartList> list, pid
// id);
inline const ParticleProperties &FindParticle(std::shared_ptr<PartList> list,
                                              pid id) {
  // position in map
  int result = -1;
  for (unsigned int i = 0; i < list->size(); ++i) {
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
                          const boost::property_tree::ptree &pt) {

  auto particleTree = pt.get_child_optional("ParticleList");
  if (!particleTree)
    return;
  for (auto const &v : particleTree.get()) {
    auto tmp = ParticleProperties(v.second);
    auto p = std::make_pair(tmp.name(), tmp);
    auto last = list->insert(p);

    if (!last.second) {
      LOG(INFO) << "ReadParticles() | Particle " << last.first->first
                << " already exists in list. We overwrite its parameters!";
      last.first->second = tmp;
    }
    tmp = last.first->second;

    // cparity is optional
    double cparity = 0.0;
    try {
      cparity = tmp.GetQuantumNumber("Cparity");
    } catch (std::exception &ex) {
    }

    LOG(DEBUG) << "ReadParticles() | Particle " << tmp.name()
               << " (id=" << tmp.GetId() << ") "
               << " J(PC)=" << tmp.GetSpinQuantumNumber("Spin") << "("
               << tmp.GetQuantumNumber("Parity") << cparity << ") "
               << " mass=" << tmp.GetMass()
               << " decayType=" << tmp.GetDecayType();
  }

  return;
}

/// Read list of particles from a boost::property_tree
void ReadParticles(PartList &list, const boost::property_tree::ptree &pt);

/// Read list of particles from a stream
inline void ReadParticles(std::shared_ptr<PartList> list,
                          std::stringstream &stream) {
  boost::property_tree::ptree tree;
  boost::property_tree::xml_parser::read_xml(stream, tree);
  ReadParticles(list, tree);
}

/// Read list of particles from a string. Note that the string contains the full
/// list; it is not a file name!
inline void ReadParticles(std::shared_ptr<PartList> list, std::string str) {
  std::stringstream ss;
  ss << str;
  ReadParticles(list, ss);
}

/// Save particle list to boost::property_tree
inline boost::property_tree::ptree
SaveParticles(std::shared_ptr<PartList> list) {
  boost::property_tree::ptree pt;

  for (auto &i : *list.get()) {
    pt.add_child("Particle", i.second.Save());
  }
  return pt;
}

/// Save particle list to boost::property_tree
inline boost::property_tree::ptree SaveParticles(PartList list) {
  boost::property_tree::ptree pt;
  for (auto &i : list) {
    pt.add_child("Particle", i.second.Save());
  }
  return pt;
}

/// Save particle list to file
inline void SaveParticles(std::shared_ptr<PartList> list,
                          std::string fileName) {
  boost::property_tree::ptree pt;
  pt.add_child("ParticleList", SaveParticles(list));
  boost::property_tree::xml_parser::write_xml(fileName, pt, std::locale());
  return;
}

inline void UpdateNode(std::shared_ptr<FitParameter> p,
                       boost::property_tree::ptree &tr) {
  for (auto &v : tr.get_child("")) {
    if (v.first == "Parameter") {
      std::string nn = v.second.get<std::string>("<xmlattr>.Name");
      std::string tt = v.second.get<std::string>("<xmlattr>.Type");
      if (nn == p->name()) {
        //        LOG(DEBUG) << "UpdateNode() | Updating node " << v.first
        //        << "." << nn << " to " << p->value();
        v.second = p->save();
        v.second.put("<xmlattr>.Type", tt);
      }
    } else {
      // Call this function recursively to find multiple matches.
      UpdateNode(p, v.second);
    }
  }
  return;
}

/// Update particle list. Note that \p partL is modified!
inline void UpdateParticleList(std::shared_ptr<PartList> &partL,
                               ParameterList &pars) {

  boost::property_tree::ptree partTr;
  partTr.add_child("ParticleList", SaveParticles(partL));
  // Loop over (double) parameters
  for (auto i : pars.doubleParameters()) {
    auto name = i->name();
    // Some default values may not have a name. We skip those.
    if (name == "")
      continue;
    // Search for parameter in tree and update it
    UpdateNode(i, partTr);
  }

  // Create new list from modified tree and assign it to \p partL
  auto newList = std::make_shared<PartList>();
  ReadParticles(newList, partTr);
  partL = newList;
  return;
}

/// Update particle list. Note that \p partL is modified!
inline void UpdateParticleList(PartList &partL, ParameterList &pars) {

  boost::property_tree::ptree partTr;
  partTr.add_child("ParticleList", SaveParticles(partL));
  // Loop over (double) parameters
  for (auto i : pars.doubleParameters()) {
    auto name = i->name();
    // Some default values may not have a name. We skip those.
    if (name == "")
      continue;
    // Search for parameter in tree and update it
    UpdateNode(i, partTr);
  }

  // Create new list from modified tree and assign it to \p partL
  PartList newList;
  ReadParticles(newList, partTr);
  partL = newList;
  return;
}

} // namespace ComPWA

#endif
