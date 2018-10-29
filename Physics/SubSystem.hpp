// Copyright (c)  2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains SubSystem class.
///

#ifndef COMPWA_PHYSICS_SUBSYSTEM_HPP_
#define COMPWA_PHYSICS_SUBSYSTEM_HPP_

#include <boost/property_tree/ptree.hpp>
#include <vector>

namespace ComPWA {
namespace Physics {

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

  SubSystem(const std::vector<std::vector<unsigned int>> &FinalStates,
            const std::vector<unsigned int> &Recoil,
            const std::vector<unsigned int> &ParentRecoil);

  virtual ~SubSystem(){};

  virtual boost::property_tree::ptree save() const;

  virtual void load(const boost::property_tree::ptree &pt);

  bool operator==(const SubSystem &b) const;

  friend std::ostream &operator<<(std::ostream &stream, const SubSystem &s);

  virtual void setFinalStates(const std::vector<std::vector<unsigned int>> &v);

  virtual const std::vector<std::vector<unsigned int>> &getFinalStates() const;

  virtual void setRecoilState(const std::vector<unsigned int> &r);

  virtual const std::vector<unsigned int> &getRecoilState() const;

  virtual void setParentRecoilState(const std::vector<unsigned int> &r);

  virtual const std::vector<unsigned int> &getParentRecoilState() const;

protected:
  std::vector<unsigned int> RecoilFinalState;
  std::vector<unsigned int> ParentRecoilFinalState;
  std::vector<std::vector<unsigned int>> FinalStates;
};

/// Helper funtions to transfor a string of space-separated numbers to a
/// vector<unsigned int>. E.g. "1 2 3" =? vector<unsigned int>({1,2,3})
inline std::vector<unsigned int> stringToVectInt(std::string str) {
  std::vector<unsigned int> result;
  std::istringstream iStr(str);
  std::vector<std::string> stringFrag{std::istream_iterator<std::string>{iStr},
                                      std::istream_iterator<std::string>{}};
  for (auto i : stringFrag) {
    result.push_back(std::stoul(i));
  }
  return result;
}

} // namespace Physics
} // namespace ComPWA

#endif
