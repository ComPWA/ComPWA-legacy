// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Properties.hpp"
#include "Core/Utils.hpp"

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
        double Mag = v.second.get<double>("<xmlattr>.Value");
        Spin J(Mag);
        try { // Projection of spin is optional (e.g. (I,I3))
          double Proj = v.second.get<double>("<xmlattr>.Projection");
          J = Spin(Mag, Proj);
        } catch (std::exception &ex) {
        }
        spinQuantumNumbers_.insert(std::make_pair(type, J));
      } else if (v.second.get<std::string>("<xmlattr>.Class") == "Int") {
        auto value = v.second.get<int>("<xmlattr>.Value");
        intQuantumNumbers_.insert(std::make_pair(type, value));
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

void ReadParticles(PartList &list, const boost::property_tree::ptree &pt) {

  auto particleTree = pt.get_child_optional("ParticleList");
  if (!particleTree)
    return;
  for (auto const &v : particleTree.get()) {
    auto tmp = ParticleProperties(v.second);
    auto p = std::make_pair(tmp.getName(), tmp);
    auto last = list.insert(p);

    if (!last.second) {
      LOG(INFO) << "ReadParticles() | Particle " << last.first->first
                << " already exists in list. We overwrite its parameters!";
      last.first->second = tmp;
    }
    tmp = last.first->second;

    // cparity is optional
    double cparity = 0.0;
    try {
      cparity = tmp.getQuantumNumber("Cparity");
    } catch (std::exception &ex) {
    }

    LOG(DEBUG) << "ReadParticles() | Particle " << tmp.getName()
               << " (id=" << tmp.getId() << ") "
               << " J(PC)=" << tmp.getSpinQuantumNumber("Spin") << "("
               << tmp.getQuantumNumber("Parity") << cparity << ") "
               << " mass=" << tmp.getMass().Value
               << " decayType=" << tmp.getDecayType();
  }

  return;
}

/*boost::property_tree::ptree ParticleProperties::save() {
  boost::property_tree::ptree pt;
  pt.put("<xmlattr>.Name", Name);
  pt.put("Pid", Id);
  pt.put("Mass", Mass);
  for (auto &i : spinQuantumNumbers_) {
    boost::property_tree::ptree tmp;
    tmp.put("<xmlattr>.Class", "Spin");
    tmp.put("<xmlattr>.Type", i.first);
    tmp.put("<xmlattr>.Value", i.second.GetSpin());
    pt.add_child("QuantumNumber", tmp);
  }
  for (auto &i : intQuantumNumbers_) {
    boost::property_tree::ptree tmp;
    tmp.put("<xmlattr>.Class", "Int");
    tmp.put("<xmlattr>.Type", i.first);
    tmp.put("<xmlattr>.Value", i.second);
    pt.add_child("QuantumNumber", tmp);
  }
  pt.add_child("DecayInfo", DecayInfo);
  return pt;
}*/

int ParticleProperties::getQuantumNumber(std::string type) const {
  auto it = intQuantumNumbers_.find(type);
  if (it == intQuantumNumbers_.end())
    throw std::runtime_error("ParticleProperties::GetQuantumNumber() | "
                             "Quantum Number '" +
                             type + "' not found!");

  return it->second;
}

ComPWA::Spin ParticleProperties::getSpinQuantumNumber(std::string type) const {
  auto it = spinQuantumNumbers_.find(type);
  if (it == spinQuantumNumbers_.end())
    throw std::runtime_error("ParticleProperties::GetSpinQuantumNumber() | "
                             "Quantum Number '" +
                             type + "' not found!");

  return it->second;
}

} // namespace ComPWA
