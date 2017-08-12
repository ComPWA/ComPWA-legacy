// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <stdio.h>
#include "Core/Properties.hpp"

using namespace ComPWA;

ParticleProperties::ParticleProperties(boost::property_tree::ptree pt)
    : Properties(pt.get<std::string>("<xmlattr>.Name"), pt.get<pid>("Pid")) {

  for (const auto &v : pt.get_child("")) {
    if (v.first == "QuantumNumber") {
      // QuantumNumbers which can be of type int or ComPWA::Spin
      std::string type = v.second.get<std::string>("<xmlattr>.Type");

      // We have to distinguish between spin and integer quantum numbers
      if (v.second.get<std::string>("<xmlattr>.Class") == "Spin") {
        auto value = v.second.get<double>("<xmlattr>.Value");
        double valueZ = 0.0;
        try { // Projection of spin is optional (e.g. (I,I3))
          valueZ = v.second.get<double>("<xmlattr>.Projection");
        } catch (std::exception &ex) {
        }
        spinQuantumNumbers_.insert(
            std::make_pair(type, ComPWA::Spin(value, valueZ)));
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
      _mass = DoubleParameterFactory(v.second);
    } else {
    }
  }

  // Info on the particle decay is stored as it is as property_tree and later
  // used by AbstractDynamicalFunctions (e.g. RelativisticBreitWigner).
  auto decayInfo = pt.get_child_optional("DecayInfo");
  if (decayInfo) {
    _decayInfo = decayInfo.get();
  } else {
    _decayInfo.put("<xmlattr>.Type", "Stable");
  }
}

int ParticleProperties::GetQuantumNumber(std::string type) const {
  auto it = intQuantumNumbers_.find(type);
  if (it == intQuantumNumbers_.end())
    throw std::runtime_error("ParticleProperties::GetQuantumNumber() | "
                             "Quantum Number '" +
                             type + "' not found!");

  return it->second;
}

ComPWA::Spin ParticleProperties::GetSpinQuantumNumber(std::string type) const {
  auto it = spinQuantumNumbers_.find(type);
  if (it == spinQuantumNumbers_.end())
    throw std::runtime_error("ParticleProperties::GetSpinQuantumNumber() | "
                             "Quantum Number '" +
                             type + "' not found!");

  return it->second;
}
