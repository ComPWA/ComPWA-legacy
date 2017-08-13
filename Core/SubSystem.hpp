// Copyright (c)  2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains SubSystem class.
///

#ifndef SubSystem_h
#define SubSystem_h

#include <vector>
#include <boost/property_tree/ptree.hpp>

namespace ComPWA {

/// \class SubSystem
/// Definition of a two-body decay node within a sequential decay tree.
/// Class contains lists for both final states of the two-body decay and a list
/// for all recoiling particles. This information is needed to calculate
/// invariant mass and angles at a two-body decay node.
class SubSystem {
public:
  SubSystem(std::vector<int> recoilS,
            std::vector<std::vector<int>> finalStates);

  /// Constructor for an isobar model in which it is assumed that particles
  /// always decay to two final states.
  SubSystem(std::vector<int> recoilS, std::vector<int> finalA,
            std::vector<int> finalB);

  virtual std::string to_string() const;

  bool operator==(const SubSystem &b) const;

  friend std::ostream &operator<<(std::ostream &stream, const SubSystem &s) {
    stream << s.to_string();
    return stream;
  }
  virtual void SetFinalStates(std::vector<std::vector<int>> v);

  virtual const std::vector<std::string> &GetFinalStatesNames() const;

  virtual void SetFinalStatesNames(std::vector<std::string> n);

  virtual const std::vector<std::vector<int>> &GetFinalStates() const;

  virtual const std::vector<int> GetHelicities() const;

  virtual void SetHelicities(std::vector<int> hel);

  virtual void SetRecoilState(const std::vector<int> r);

  virtual const std::vector<int> &GetRecoilState() const;

protected:
  std::string _title;
  std::vector<int> _recoilState;
  std::vector<std::vector<int>> _finalStates;
  std::vector<int> helicities_;
  std::vector<std::string> _finalStatesNames;
};

/// Helper funtions to transfor a string of space-separated numbers to a
/// vector<int>. E.g. "1 2 3" =? vector<int>({1,2,3})
inline std::vector<int> stringToVectInt(std::string str) {
  std::vector<int> result;
  std::istringstream iStr(str);
  std::vector<std::string> stringFrag{std::istream_iterator<std::string>{iStr},
                                      std::istream_iterator<std::string>{}};
  for (auto i : stringFrag) {
    result.push_back(std::stoi(i));
  }
  return result;
}

inline SubSystem SubSystemFactory(const boost::property_tree::ptree pt) {
  // Read subSystem definition
  std::vector<int> recoilState;
  auto recoil =
      pt.get_optional<std::string>("RecoilSystem.<xmlattr>.FinalState");
  if (recoil) {
    recoilState = stringToVectInt(recoil.get());
  }

  auto decayProducts = pt.get_child("DecayProducts");
  std::vector<std::vector<int>> finalStates;
  std::vector<int> helicities;
  std::vector<std::string> finalStatesNames;
  for (auto i : decayProducts) {
    finalStates.push_back(
        stringToVectInt(i.second.get<std::string>("<xmlattr>.FinalState")));
    finalStatesNames.push_back(i.second.get<std::string>("<xmlattr>.Name"));
    helicities.push_back(i.second.get<int>("<xmlattr>.Helicity"));
  }
  SubSystem subSys(recoilState, finalStates);
  subSys.SetFinalStatesNames(finalStatesNames);
  subSys.SetHelicities(helicities);
  return subSys;
}

inline boost::property_tree::ptree SubSystemSave(SubSystem sys) {
  boost::property_tree::ptree pt;
  auto recoilV = sys.GetRecoilState();
  if (recoilV.size()) {
    std::string recoilStr;
    for (auto i = 0; i < recoilV.size(); i++) {
      recoilStr += std::to_string(recoilV.at(i));
      if (i < recoilV.size() - 1)
        recoilStr += " ";
    }
    pt.put("RecoilSystem.<xmlattr>.FinalState", recoilStr);
  }

  boost::property_tree::ptree daughterTr;

  // Information daugher final state A
  auto finalS = sys.GetFinalStates();
  auto helicities = sys.GetHelicities();
  auto finalSNames = sys.GetFinalStatesNames();
  for (int j = 0; j < finalS.size(); j++) {
    std::string strA;
    if (finalS.at(j).size()) {
      for (auto i = 0; i < finalS.at(j).size(); i++) {
        strA += std::to_string(finalS.at(j).at(i));
        if (i < finalS.at(j).size() - 1)
          strA += " ";
      }
    }

    boost::property_tree::ptree fTr;
    fTr.put("<xmlattr>.Name", finalSNames.at(j));
    fTr.put("<xmlattr>.FinalState", strA);
    fTr.put("<xmlattr>.Helicity", helicities.at(j));
    daughterTr.add_child("Particle", fTr);
  }

  pt.add_child("DecayProducts", daughterTr);

  return pt;
}

} /* namespace ComPWA */

#endif /* SubSystem_h */
