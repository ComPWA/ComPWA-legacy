// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PROPERTIES_HPP_
#define COMPWA_PROPERTIES_HPP_

#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/Exceptions.hpp"
#include "Core/FitParameter.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {

/// Particle ID.
/// Usually the pid's from PDG are used here:
/// http://pdg.lbl.gov/mc_particleId_contents.html
typedef int pid;

///
/// \class ParticleProperties
///
class ParticleProperties {
public:
  ParticleProperties(boost::property_tree::ptree pt);

  std::string getName() const { return Name; }

  pid getId() const { return Id; }

  ComPWA::FitParameter<double> getMass() const { return Mass; }

  template <typename T> T getQuantumNumber(std::string type) const;

  boost::property_tree::ptree getDecayInfo() const { return DecayInfo; }

  std::string getDecayType() const {
    return DecayInfo.get<std::string>("<xmlattr>.Type");
  }

private:
  std::string Name;
  pid Id;
  ComPWA::FitParameter<double> Mass;
  std::map<std::string, int> IntQuantumNumbers;
  std::map<std::string, double> RealQuantumNumbers;

  /// Store decay info in property_tree. The tree is later on passed to the
  /// respective class.
  boost::property_tree::ptree DecayInfo;
};

template <>
inline int ParticleProperties::getQuantumNumber(std::string type) const {
  auto it = IntQuantumNumbers.find(type);
  if (it == IntQuantumNumbers.end())
    throw std::runtime_error("ParticleProperties::getQuantumNumber<int>() | "
                             "Quantum Number '" +
                             type + "' not found!");

  return it->second;
}

template <>
inline double ParticleProperties::getQuantumNumber(std::string type) const {
  auto it = RealQuantumNumbers.find(type);
  if (it == RealQuantumNumbers.end())
    throw std::runtime_error("ParticleProperties::getQuantumNumber<double>() | "
                             "Quantum Number '" +
                             type + "' not found!");

  return it->second;
}

/// A map of particle properties is used everywhere where particle information
/// is needed. Properties are accessed by the particle name.
/// Note: Propably would be better to access particles by their pid?
typedef std::map<std::string, ParticleProperties> PartList;

inline std::ostream &operator<<(std::ostream &os, const PartList &p) {
  for (auto i : p)
    os << i.first << " [ " << i.second.getId()
       << " ]: mass = " << i.second.getMass().Value << std::endl;
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
    if (it->second.getId() == id) {
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
    auto p = std::make_pair(tmp.getName(), tmp);
    auto last = list->insert(p);

    if (!last.second) {
      LOG(INFO) << "ReadParticles() | Particle " << last.first->first
                << " already exists in list. We overwrite its parameters!";
      last.first->second = tmp;
    }
    tmp = last.first->second;

    // cparity is optional
    int cparity = 0;
    try {
      cparity = tmp.getQuantumNumber<int>("Cparity");
    } catch (std::exception &ex) {
    }

    LOG(DEBUG) << "ReadParticles() | Particle " << tmp.getName()
               << " (id=" << tmp.getId() << ") "
               << " J(PC)=" << tmp.getQuantumNumber<double>("Spin") << "("
               << tmp.getQuantumNumber<int>("Parity") << cparity << ") "
               << " mass=" << tmp.getMass().Value
               << " decayType=" << tmp.getDecayType();
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

} // namespace ComPWA

#endif
