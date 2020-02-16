// Copyright (c)  2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains SubSystem class.
///

#ifndef COMPWA_PHYSICS_SUBSYSTEM_HPP_
#define COMPWA_PHYSICS_SUBSYSTEM_HPP_

#include <functional>
#include <vector>

#include <boost/property_tree/ptree.hpp>

namespace ComPWA {

typedef std::vector<unsigned int> IndexList;

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
  SubSystem(const std::vector<std::vector<unsigned int>> &FinalStates,
            const std::vector<unsigned int> &Recoil,
            const std::vector<unsigned int> &ParentRecoil);

  virtual ~SubSystem() = default;

  bool operator==(const SubSystem &b) const;

  friend std::ostream &operator<<(std::ostream &stream, const SubSystem &s);
  virtual const std::vector<std::vector<unsigned int>> &getFinalStates() const;
  virtual const std::vector<unsigned int> &getRecoilState() const;
  virtual const std::vector<unsigned int> &getParentRecoilState() const;

private:
  std::vector<std::vector<unsigned int>> DecayProductFinalStates;
  std::vector<unsigned int> RecoilFinalState;
  std::vector<unsigned int> ParentRecoilFinalState;
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

namespace std {

// this function was copied from boost
template <typename T> std::size_t combineHash(size_t hash, const T &v) {
  hash ^= std::hash<T>()(v) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
  return hash;
}
template <> struct hash<ComPWA::IndexList> {
  std::size_t operator()(const ComPWA::IndexList &key) const {
    std::size_t seed = key.size();
    for (unsigned int i : key) {
      seed = combineHash(seed, i);
    }
    return seed;
  }
};

template <> struct hash<ComPWA::Physics::SubSystem> {
  std::size_t operator()(const ComPWA::Physics::SubSystem &key) const {
    std::size_t seed = key.getFinalStates().size();
    for (auto const &x : key.getFinalStates()) {
      seed = combineHash(seed, x);
    }
    seed = combineHash(seed, key.getRecoilState());
    seed = combineHash(seed, key.getParentRecoilState());

    return seed;
  }
};

} // namespace std

#endif
