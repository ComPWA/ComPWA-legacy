// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Properties.hpp"

#include "Core/Logging.hpp"
#include "Core/Utils.hpp"

#include "boost/property_tree/xml_parser.hpp"

namespace ComPWA {

ParticleProperties::ParticleProperties(boost::property_tree::ptree pt) {
  Name = pt.get<std::string>("<xmlattr>.Name");
  Id = pt.get<pid>("Pid");
  for (const auto &v : pt.get_child("")) {
    if (v.first == "QuantumNumber") {
      // QuantumNumbers which can be of type int or ComPWA::Spin
      std::string type = v.second.get<std::string>("<xmlattr>.Type");

      // We have to distinguish between spin and integer quantum numbers
      if (v.second.get<std::string>("<xmlattr>.Class") == "Spin") {
        auto value = v.second.get<double>("<xmlattr>.Value");
        try { // Projection of spin is optional (e.g. (I,I3))
          double valueZ = v.second.get<double>("<xmlattr>.Projection");
          RealQuantumNumbers.insert(
              std::make_pair(type + ".Projection", valueZ));
        } catch (std::exception &ex) {
        }
        RealQuantumNumbers.insert(std::make_pair(type, value));
      } else if (v.second.get<std::string>("<xmlattr>.Class") == "Int") {
        auto value = v.second.get<int>("<xmlattr>.Value");
        IntQuantumNumbers.insert(std::make_pair(type, value));
      } else {
        throw BadParameter(
            "ParticleProperties::ParticleProperties() | "
            "QuantumNumber is neither of type 'Spin' nor of type "
            "'Int'!");
      }
    } else if (v.first == "Parameter") {
      // Parameter (e.g. Mass)
      if (v.second.get<std::string>("<xmlattr>.Type") != "Mass")
        continue;
      Mass.Name = v.second.get<std::string>("<xmlattr>.Name");
      Mass.Value = v.second.get<double>("Value");

      auto error = v.second.get_optional<double>("Error");
      if (error) {
        Mass.Error = std::make_pair(error.get(), error.get());
      } else {
        auto errorlow = v.second.get_optional<double>("ErrorLow");
        auto errorhigh = v.second.get_optional<double>("ErrorHigh");
        if (errorlow && errorhigh) {
          Mass.Error = std::make_pair(errorlow.get(), errorhigh.get());
        } else if (errorlow) {
          throw std::runtime_error("ParticleProperties | Error of Parameter "
                                   "Mass not properly set (missing ErrorHigh)");
        } else if (errorhigh) {
          throw std::runtime_error("ParticleProperties | Error of Parameter "
                                   "Mass not properly set (missing ErrorLow)");
        }
      }

      auto fix = v.second.get_optional<bool>("Fix");
      if (fix)
        Mass.IsFixed = fix.get();
      else
        Mass.IsFixed = true;

      auto min = v.second.get_optional<double>("Min");
      auto max = v.second.get_optional<double>("Max");
      if (min && max) {
        Mass.HasBounds = true;
        Mass.Bounds = std::make_pair(min.get(), max.get());
      } else if (min || max) { // Bounds not completely specified
        throw std::runtime_error(
            "FitParameterFactory() | Parameter bounds not properly set!");
      }
    } else {
    }
  }

  // Info on the particle decay is stored as it is as property_tree and later
  // used by AbstractDynamicalFunctions (e.g. RelativisticBreitWigner).
  auto decayInfo = pt.get_child_optional("DecayInfo");
  if (decayInfo) {
    DecayInfo = decayInfo.get();
  } else {
    DecayInfo.put("<xmlattr>.Type", "Stable");
  }
}

/// Read list of particles from a boost::property_tree
void insertParticles(ParticleList &list,
                     const boost::property_tree::ptree &pt) {
  auto particleTree = pt.get_child_optional("ParticleList");
  if (!particleTree)
    return;
  for (auto const &v : particleTree.get()) {
    if (v.first != "Particle")
      continue;
    ParticleProperties tmp(v.second);
    auto result = list.insert(tmp);

    if (!result.second) {
      LOG(INFO) << " Particle " << tmp.getName() << " with identical ID "
                << tmp.getId() << " already exists in list with the name "
                << result.first->getName() << " and ID "
                << result.first->getId()
                << ". Particle properties will be overwritten!";
      list.erase(result.first);
      result = list.insert(tmp);
    }
    tmp = *result.first;

    std::stringstream ss;
    ss << " Particle " << tmp.getName() << " (id=" << tmp.getId()
       << ") J(PC)=" << tmp.getQuantumNumber<double>("Spin") << "(";
    // parities are optional
    try {
      int parity = tmp.getQuantumNumber<int>("Parity");
      ss << parity;
    } catch (std::exception &ex) {
    }
    try {
      int cparity = tmp.getQuantumNumber<int>("Cparity");
      ss << cparity;
    } catch (std::exception &ex) {
    }
    ss << ") mass=" << tmp.getMass().Value
       << " decayType=" << tmp.getDecayType();
    LOG(DEBUG) << ss.str();
  }

  return;
}

void insertParticles(ParticleList &list, std::stringstream &Stream) {
  boost::property_tree::ptree tree;
  boost::property_tree::xml_parser::read_xml(Stream, tree);
  insertParticles(list, tree);
}

void insertParticles(ParticleList &list, std::string FileName) {
  boost::property_tree::ptree tree;
  boost::property_tree::xml_parser::read_xml(FileName, tree);
  insertParticles(list, tree);
}

ParticleList readParticles(std::stringstream &Stream) {
  boost::property_tree::ptree tree;
  boost::property_tree::xml_parser::read_xml(Stream, tree);
  return readParticles(tree);
}

ParticleList readParticles(std::string FileName) {
  boost::property_tree::ptree tree;
  boost::property_tree::xml_parser::read_xml(FileName, tree);
  return readParticles(tree);
}

ParticleList readParticles(boost::property_tree::ptree &tree) {
  ParticleList list;
  insertParticles(list, tree);
  return list;
}

} // namespace ComPWA
