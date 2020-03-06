// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PROPERTIES_HPP_
#define COMPWA_PROPERTIES_HPP_

#include "Core/Exceptions.hpp"
#include "Core/FitParameter.hpp"

#include "boost/property_tree/ptree.hpp"

#include <map>
#include <set>
#include <vector>

namespace ComPWA {

/// Particle ID.
/// Usually the PIDs from PDG are used here:
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

  friend bool operator<(const ParticleProperties &l,
                        const ParticleProperties &r) {
    return l.Id < r.Id;
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

using ParticleList = std::set<ParticleProperties>;

inline std::ostream &operator<<(std::ostream &os, const ParticleList &p) {
  for (auto i : p)
    os << i.getName() << " [ " << i.getId()
       << " ]: mass = " << i.getMass().Value << std::endl;
  return os;
}

inline const ParticleProperties &findParticle(const ParticleList &list,
                                              pid Pid) {
  auto found = std::find_if(list.begin(), list.end(),
                            [&Pid](auto const &x) { return x.getId() == Pid; });
  if (list.end() == found) {
    throw std::runtime_error("Could not find particle with id " +
                             std::to_string(Pid) + " in list");
  }
  return *found;
}

inline const ParticleProperties &findParticle(const ParticleList &list,
                                              std::string refname) {
  auto found =
      std::find_if(list.begin(), list.end(), [&refname](auto const &x) {
        return x.getName() == refname;
      });
  if (list.end() == found) {
    throw std::runtime_error("Could not find particle with name " + refname +
                             " in list");
  }
  return *found;
}

/// insert particles from a boost::property_tree into a ParticleList
void insertParticles(ParticleList &list, const boost::property_tree::ptree &pt);

/// insert particles from a stringstream into a ParticleList
void insertParticles(ParticleList &list, std::stringstream &Stream);

/// insert particles from a xml file into a ParticleList
void insertParticles(ParticleList &list, std::string FileName);

/// Read list of particles from a stringstream
/// For some reason the boost xml parser needs a non-const reference
ParticleList readParticles(std::stringstream &Stream);

/// Read list of particles from a xml file
ParticleList readParticles(std::string FileName);

/// Read list of particles from a `boost::property_tree::ptree`
ParticleList readParticles(boost::property_tree::ptree &pt);

} // namespace ComPWA

#endif
