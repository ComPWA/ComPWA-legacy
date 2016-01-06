//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef PHYSICS_HELICITYAMPLITUDE_PARTICLESTATEDEFINITIONS_HPP_
#define PHYSICS_HELICITYAMPLITUDE_PARTICLESTATEDEFINITIONS_HPP_

#include <set>
#include <map>

#include "boost/property_tree/ptree.hpp"

#include "Core/Utility.hpp"

namespace HelicityFormalism {

typedef boost::property_tree::ptree DynamicalInfo;

struct HelicityAngles {
  double theta_;
  double phi_;

  HelicityAngles(double theta, double phi) :
      theta_(theta), phi_(phi) {
  }
};

struct IDInfo {
  int particle_id_;
  std::string name_;

  bool operator==(const IDInfo &rhs) const {
    /* if (this->id_ != rhs.id_)
     return false;*/
    if (this->particle_id_ != rhs.particle_id_)
      return false;
    if (this->name_ != rhs.name_)
      return false;

    return true;
  }
  bool operator!=(const IDInfo &rhs) const {
    return !(*this == rhs);
  }

  bool operator<(const IDInfo &rhs) const {
    /*if (this->id_ < rhs.id_)
     return true;
     else if (this->id_ > rhs.id_)
     return false;*/

    return lessThenIgnoringID(*this, rhs);
  }

  static bool lessThenIgnoringID(const IDInfo &lhs, const IDInfo &rhs) {
    if (lhs.particle_id_ < rhs.particle_id_)
      return true;
    else if (lhs.particle_id_ > rhs.particle_id_)
      return false;
    if (lhs.name_ < rhs.name_)
      return true;

    return false;
  }

  bool operator>(const IDInfo &rhs) const {
    return (rhs < *this);
  }
};



struct ParticleStateInfo {
  unsigned int unique_id_;
  IDInfo pid_information_;
  ComPWA::Spin spin_information_;
  DynamicalInfo dynamical_information_;

  bool operator==(const ParticleStateInfo &rhs) const {
    if (this->unique_id_ != rhs.unique_id_)
      return false;
    if (this->pid_information_ != rhs.pid_information_)
      return false;
    if (this->spin_information_ != rhs.spin_information_)
      return false;

    return true;
  }

  bool operator!=(const ParticleStateInfo &rhs) const {
    return !((*this) == rhs);
  }

  bool operator<(const ParticleStateInfo &rhs) const {
    if (this->unique_id_ < rhs.unique_id_)
      return true;
    else if (this->unique_id_ > rhs.unique_id_)
      return false;
    if (this->pid_information_ < rhs.pid_information_)
      return true;
    else if (this->pid_information_ > rhs.pid_information_)
      return false;
    if (this->spin_information_ < rhs.spin_information_)
      return true;

    return false;
  }
  bool operator>(const ParticleStateInfo &rhs) const {
    return (rhs < *this);
  }
};

struct ParticleIDComparison {
  unsigned int ps_id_;

  ParticleIDComparison() {
  }
  ParticleIDComparison(unsigned int ps_id) :
      ps_id_(ps_id) {
  }

  bool operator()(const ParticleStateInfo& ps) {
    return ps.unique_id_ == ps_id_;
  }

  bool operator()(const ParticleStateInfo& lhs, const ParticleStateInfo& rhs) {
    return lhs.unique_id_ < rhs.unique_id_;
  }
};

typedef std::pair<ComPWA::Spin, ComPWA::Spin> SpinInfoPair;
typedef std::pair<DynamicalInfo, DynamicalInfo> DynamicalInfoPair;

typedef std::vector<unsigned int> IndexList;
typedef std::pair<unsigned int, unsigned int> IndexPair;
typedef std::map<unsigned int, unsigned int> IndexMapping;

struct TwoBodyDecayIndices {
  unsigned int mother_index_;
  unsigned int decay_state_index_;
  IndexPair decay_products_;

  bool operator<(const TwoBodyDecayIndices &rhs) const {
    if (this->mother_index_ < rhs.mother_index_)
      return true;
    else if (this->mother_index_ > rhs.mother_index_)
      return false;
    if (this->decay_state_index_ < rhs.decay_state_index_)
      return true;
    else if (this->decay_state_index_ > rhs.decay_state_index_)
      return false;
    if (this->decay_products_ < rhs.decay_products_)
      return true;

    return false;
  }

  bool operator==(const TwoBodyDecayIndices &rhs) const {
    if (this->mother_index_ != rhs.mother_index_)
      return false;
    if (this->decay_state_index_ != rhs.decay_state_index_)
      return false;
    if (this->decay_products_ != rhs.decay_products_)
      return false;

    return true;
  }

};

struct TwoBodyDecaySpinInformation {
  ComPWA::Spin initial_state_;
  SpinInfoPair final_state_;

  bool operator<(const TwoBodyDecaySpinInformation &rhs) const {
    if (this->initial_state_ < rhs.initial_state_)
      return true;
    else if (this->initial_state_ > rhs.initial_state_)
      return false;
    if (this->final_state_ < rhs.final_state_)
      return true;

    return false;
  }

  bool operator>(const TwoBodyDecaySpinInformation &rhs) const {
    return (rhs < *this);
  }
};

struct TwoBodyDecayDynamicalInformation {
  DynamicalInfo initial_state_;

  bool operator<(const TwoBodyDecayDynamicalInformation &rhs) const {
    if (generateMap(this->initial_state_) < generateMap(rhs.initial_state_))
      return true;

    return false;
  }

  bool operator>(const TwoBodyDecayDynamicalInformation &rhs) const {
    return (rhs < *this);
  }

  std::map<std::string, std::string> generateMap(
      const DynamicalInfo& tree) const {
    std::map<std::string, std::string> return_map;
    DynamicalInfo::const_iterator tree_iter;
    for (tree_iter = tree.begin(); tree_iter != tree.end(); ++tree_iter) {
      return_map[tree_iter->first] = tree_iter->second.get_value<std::string>();
    }
    return return_map;
  }
};

struct TwoBodyDecayInformation {
  TwoBodyDecaySpinInformation spin_info_;
  TwoBodyDecayDynamicalInformation dynamical_info_;

  bool operator<(const TwoBodyDecayInformation &rhs) const {
    if (this->spin_info_ < rhs.spin_info_)
      return true;
    else if (this->spin_info_ > rhs.spin_info_)
      return false;
    if (this->dynamical_info_ < rhs.dynamical_info_)
      return true;

    return true;
  }
};

struct checkLessThanOnIDInfoVectorsIgnoringID {
  bool operator()(const std::vector<IDInfo>& lhs,
      const std::vector<IDInfo>& rhs) const {
    if (lhs.size() < rhs.size())
      return true;
    else if (lhs.size() > rhs.size())
      return false;
    std::vector<unsigned int> not_found_elements;
    std::vector<IDInfo> clone_lhs(lhs);
    std::vector<IDInfo> clone_rhs(rhs);
    std::sort(clone_lhs.begin(), clone_lhs.end(), IDInfo::lessThenIgnoringID);
    std::sort(clone_rhs.begin(), clone_rhs.end(), IDInfo::lessThenIgnoringID);
    for (unsigned int i = 0; i < lhs.size(); ++i) {
      if (clone_lhs[i].particle_id_ < clone_rhs[i].particle_id_)
        return true;
      else if (clone_lhs[i].particle_id_ > clone_rhs[i].particle_id_)
        return false;
    }
    return false;
  }
};

struct TwoBodyDecayTopology {
  IDInfo top_node_id_info_;
  std::vector<IDInfo> final_state_id_list_;
  // list of final state particle lists corresponding to decay states
  // this actually uniquely defines a decay topology
  // all the other information is stored for other purposes
  std::vector<std::vector<IDInfo> > final_state_content_id_lists_;
  std::vector<TwoBodyDecayIndices> decay_node_fs_content_index_infos_;

  // evaluation order vector
  IndexList unique_id_decay_node_order_;

  // this stuff is for final state combinatorics
  unsigned int top_node_unique_id_;
  std::map<unsigned int, IndexList> particle_unique_id_decay_tree_;
  std::map<unsigned int, IndexList> final_state_content_unique_id_mapping_;

  unsigned int insertFinalStateContentList(
      std::vector<IDInfo> fs_content_list) {
    std::sort(fs_content_list.begin(), fs_content_list.end());

    unsigned int position_index;
    auto result = std::find(final_state_content_id_lists_.begin(),
        final_state_content_id_lists_.end(), fs_content_list);
    if (result != final_state_content_id_lists_.end())
      position_index = result - final_state_content_id_lists_.begin();
    else {
      position_index = final_state_content_id_lists_.size();
      final_state_content_id_lists_.push_back(fs_content_list);
    }

    return position_index;
  }

  unsigned int findMotherIndex(unsigned int decay_state_index) const {
    unsigned int mother_index(decay_state_index);
    for (auto iter = particle_unique_id_decay_tree_.begin();
        iter != particle_unique_id_decay_tree_.end(); ++iter) {
      if (std::find(iter->second.begin(), iter->second.end(), decay_state_index)
          != iter->second.end()) {
        mother_index = iter->first;
        break;
      }
    }
    return mother_index;
  }

  unsigned int findFinalStateList(std::vector<IDInfo> fs_content_list) {
    std::sort(fs_content_list.begin(), fs_content_list.end());

    unsigned int position_index(0);
    auto result = std::find(final_state_content_id_lists_.begin(),
        final_state_content_id_lists_.end(), fs_content_list);
    if (result != final_state_content_id_lists_.end())
      position_index = result - final_state_content_id_lists_.begin();

    return position_index;
  }

  void setMotherDecayIndexForDecayIndex(unsigned int mother_index,
      unsigned int decay_index) {
    auto result =
        std::find_if(decay_node_fs_content_index_infos_.begin(),
            decay_node_fs_content_index_infos_.end(),
            [&](const TwoBodyDecayIndices tbdi) {return tbdi.decay_state_index_ == decay_index;});
    if (result != decay_node_fs_content_index_infos_.end())
      result->mother_index_ = mother_index;
  }

  bool operator==(const TwoBodyDecayTopology &rhs) const {
    std::vector<std::vector<IDInfo> > fscl_copy_lhs(
        this->final_state_content_id_lists_);
    std::vector<std::vector<IDInfo> > fscl_copy_rhs(
        rhs.final_state_content_id_lists_);
    std::sort(fscl_copy_lhs.begin(), fscl_copy_lhs.end());
    std::sort(fscl_copy_rhs.begin(), fscl_copy_rhs.end());

    if (fscl_copy_lhs != fscl_copy_rhs)
      return false;

    return true;
  }

  bool operator<(const TwoBodyDecayTopology &rhs) const {
    std::vector<std::vector<IDInfo> > fscl_copy_lhs(
        this->final_state_content_id_lists_);
    std::vector<std::vector<IDInfo> > fscl_copy_rhs(
        rhs.final_state_content_id_lists_);

    std::sort(fscl_copy_lhs.begin(), fscl_copy_lhs.end());
    std::sort(fscl_copy_rhs.begin(), fscl_copy_rhs.end());

    if (fscl_copy_lhs < fscl_copy_rhs)
      return true;

    return false;
  }

  void print() const {
    /*for (auto iter = decay_node_fs_content_index_infos_.begin();
     iter != decay_node_fs_content_index_infos_.end(); ++iter) {
     std::cout << iter->decay_state_index_ << " "
     << getFSParticleListString(iter->decay_state_index_) << std::endl;
     std::cout << "originates from " << iter->mother_index_ << " "
     << getFSParticleListString(iter->mother_index_) << std::endl;
     std::cout << "decays into " << iter->decay_products_.first << " "
     << getFSParticleListString(iter->decay_products_.first) << " + "
     << iter->decay_products_.second << " "
     << getFSParticleListString(iter->decay_products_.second) << std::endl;
     }*/

    for (unsigned int i = 0; i < final_state_content_id_lists_.size(); ++i) {
      std::cout << getFSParticleListString(i) << std::endl;
    }
    std::cout << std::endl;
  }

  std::string getFSParticleListString(unsigned int index) const {
    std::stringstream ss;
    for (unsigned int i = 0; i < final_state_content_id_lists_[index].size();
        ++i) {
      ss << final_state_content_id_lists_[index][i].name_ << " ";
    }
    return ss.str();
  }
}
;

}

#endif /* PHYSICS_HELICITYAMPLITUDE_PARTICLESTATEDEFINITIONS_HPP_ */
