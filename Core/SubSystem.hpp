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

///
/// \class SubSystem
/// Definition of a two-body decay node within a sequential decay tree.
/// Class contains lists for both final states of the two-body decay and a list
/// for all recoiling particles. This information is needed to calculate
/// invariant mass and angles at a two-body decay node.
///
class SubSystem {
public:
  SubSystem(const boost::property_tree::ptree &pt);
  
  SubSystem(std::vector<int> recoilS,
            std::vector<std::vector<int>> finalStates);

  /// Constructor for an isobar model in which it is assumed that particles
  /// always decay to two final states.
  SubSystem(std::vector<int> recoilS, std::vector<int> finalA,
            std::vector<int> finalB);

  virtual boost::property_tree::ptree save() const;
  
  virtual void load(const boost::property_tree::ptree &pt);

  // Create a unique title
  virtual std::string to_string() const;

  bool operator==(const SubSystem &b) const;

  friend std::ostream &operator<<(std::ostream &stream, const SubSystem &s) {
    return stream << s.to_string();
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

} // ns::ComPWA

#endif
