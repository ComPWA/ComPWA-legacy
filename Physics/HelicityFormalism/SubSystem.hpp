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
    for (auto j=_finalStates.begin(); j != _finalStates.end(); ++j) {
      for (auto i : *j)
        stream << std::to_string(i);
      if( j != _finalStates.end()-1 ) stream << ")+(";
      else stream << ")";
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

  virtual void SetRecoilState(const std::vector<int> r) { _recoilState = r; }

  virtual const std::vector<int> &GetRecoilState() const {
    return _recoilState;
  }

protected:
  std::string title;
  std::vector<int> _recoilState;
  std::vector<std::vector<int>> _finalStates;
  std::vector<std::string> _finalStatesNames;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* SubSystem_h */
