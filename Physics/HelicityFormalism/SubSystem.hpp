//
//  SubSystem.hpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 03.08.17.
//
//

#ifndef SubSystem_h
#define SubSystem_h

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

/*! Definition of a two-body decay node within a sequential decay tree.
 Class contains lists for both final states of the two-body decay and a list
 for all recoiling particles. This information is needed to calculate
 invariant mass and angles at a two-body decay node.
 */
class SubSystem {
public:
  SubSystem(){};

  SubSystem(std::vector<int> recoilS, std::vector<std::vector<int>> finalStates)
      : _recoilState(recoilS), _finalStates(finalStates) {

    title = to_string();
    // LOG(trace) << "SubSystem::SubSystem() | Creating sub system "<<title;
  }

  SubSystem(std::vector<int> recoilS, std::vector<int> finalA,
            std::vector<int> finalB)
      : _recoilState(recoilS) {
    std::vector<std::vector<int>> tmp;
    tmp.push_back(finalA);
    tmp.push_back(finalB);
    SetFinalStates(tmp);

    title = to_string();
    // LOG(trace) << "SubSystem::SubSystem() | Creating sub system "<<title;
  }

  virtual std::string to_string() const {
    // Creating unique title
    std::stringstream stream;
    stream << "(";
    for (auto i : _recoilState)
      stream << std::to_string(i);
    stream << ")->(";
    for (auto j = _finalStates.begin(); j != _finalStates.end(); ++j) {
      for (auto i : *j)
        stream << std::to_string(i);
      if (j != _finalStates.end() - 1)
        stream << ")+(";
      else
        stream << ")";
    }

    return stream.str();
  }

  bool operator==(const SubSystem &b) const {
    if (_recoilState == b._recoilState && _finalStates == b._finalStates)
      return true;
    return false;
  }

  friend std::ostream &operator<<(std::ostream &stream, const SubSystem &s) {
    stream << s.to_string();
    return stream;
  }

  virtual void SetFinalStates(std::vector<std::vector<int>> v) {
    _finalStates = v;
  }

  virtual const std::vector<std::string> &GetFinalStatesNames() const {
    return _finalStatesNames;
  }

  virtual void SetFinalStatesNames(std::vector<std::string> n) {
    if (n.size() != _finalStates.size()) {
      throw std::runtime_error("SubSystem::SetFinalStatesNames() | Length of "
                               "vectors does not match with the number of "
                               "final states.");
    }
    _finalStatesNames = n;
  }

  virtual const std::vector<std::vector<int>> &GetFinalStates() const {
    return _finalStates;
  }

  virtual const std::vector<int> GetHelicities() const {
    if (helicities_.size() != _finalStates.size())
      throw std::runtime_error("SubSystem::GetHelicities() | Helicities are "
                               "not defined for all final states!");
    return helicities_;
  }
  virtual void SetHelicities(std::vector<int> hel) {
    if (hel.size() != _finalStates.size())
      throw std::runtime_error("SubSystem::SetHelicities() | Helicities are "
                               "not defined for all final states!");
    helicities_ = hel;
  }

  virtual void SetRecoilState(const std::vector<int> r) { _recoilState = r; }

  virtual const std::vector<int> &GetRecoilState() const {
    return _recoilState;
  }

protected:
  std::string title;
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

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* SubSystem_h */
