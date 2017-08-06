//
//  Properties.hpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 06.08.17.
//
//

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

  ParticleProperties(boost::property_tree::ptree pt)
      : Properties(pt.get<std::string>("<xmlattr>.Name"), pt.get<pid>("Pid")) {

    _mass = DoubleParameterFactory(pt.get_child("Mass"));
    // Parameters optional. Shall we require them?

    // Loop over QuantumNumbers which can be of type int or ComPWA::Spin
    for (const auto &v : pt.get_child("")) {
      if (v.first != "QuantumNumber")
        continue;
      std::string nn = v.second.get<std::string>("<xmlattr>.Name");

      // We have to distinguish between spin and integer quantum numbers
      if (v.second.get<std::string>("<xmlattr>.Type") == "Spin") {
        auto value = v.second.get<double>("<xmlattr>.Value");
        auto valueZ = v.second.get<double>("<xmlattr>.Projection");
        spinQuantumNumbers_.insert(
            std::make_pair(nn, ComPWA::Spin(value, valueZ)));
      } else if (v.second.get<std::string>("<xmlattr>.Type") == "Spin") {
        auto value = v.second.get<int>("<xmlattr>.Value");
        intQuantumNumbers_.insert(std::make_pair(nn, value));
      } else {
        throw BadParameter("ParticleProperties::ParticleProperties() | "
                          "QuantumNumber is neither of type 'Spin' nor of type "
                          "'Int'!");
      }
    }
    SetSpin(ComPWA::Spin(pt.get<double>("Spin", 0.)));
    SetIsoSpin(ComPWA::Spin(pt.get<double>("IsoSpin", 0.)));
    SetIsoSpinZ(ComPWA::Spin(pt.get<double>("IsoSpinZ", 0.)));
    SetParity(pt.get<double>("Parity", 0));
    SetCparity(pt.get<double>("Cparity", 0));
    SetGparity(pt.get<double>("Gparity", 0));

    auto decayInfo = pt.get_child_optional("DecayInfo");
    if (decayInfo)
      _decayInfo = decayInfo.get();
    else {
      _decayInfo.put("<xmlattr>.Type", "stable");
    }
  }

  double GetMass() const { return _mass.GetValue(); }

  ComPWA::DoubleParameter GetMassPar() const { return _mass; }

  int GetQuantumNumber(std::string type) const {
    auto it = intQuantumNumbers_.find(type);
    if (it == intQuantumNumbers_.end())
      throw std::runtime_error("ParticleProperties::GetQuantumNumber() | "
                               "Quantum Number not found!");

    return it->second;
  }

  ComPWA::Spin GetSpinQuantumNumber(std::string type) const {
    auto it = spinQuantumNumbers_.find(type);
    if (it == spinQuantumNumbers_.end())
      throw std::runtime_error("ParticleProperties::GetSpinQuantumNumber() | "
                               "Quantum Number not found!");

    return it->second;
  }

  boost::property_tree::ptree GetDecayInfo() const { return _decayInfo; }

  std::string GetDecayType() const {
    return _decayInfo.get<std::string>("<xmlattr>.Type");
  }
  
  //Obsolete
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
  std::map<std::string, int> intQuantumNumbers_;
  std::map<std::string, ComPWA::Spin> spinQuantumNumbers_;

  ComPWA::Spin _spin;
  ComPWA::Spin _isoSpin;
  ComPWA::Spin _isoSpinZ;
  int _parity;
  int _cparity;
  int _gparity;

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
