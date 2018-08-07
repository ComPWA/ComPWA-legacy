// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//
// \file
// Kinematics interface class
//

#ifndef KINEMATICS_HPP_
#define KINEMATICS_HPP_

#include <string>
#include <vector>

#include "Core/Particle.hpp"
#include "Core/Properties.hpp"

namespace ComPWA {

class DataPoint;
class Event;
class SubSystem;

struct KinematicsProperties {
  std::vector<pid> InitialState;
  std::vector<pid> FinalState;
  std::shared_ptr<PartList> ParticleList;
  /// Four momentum of the initial particle reaction
  ComPWA::FourMomentum InitialStateP4;
  // we use a vector instead of a map here, due to cache optimizations
  std::vector<unsigned int> FinalStateEventPositionMapping;

  KinematicsProperties() {}

  KinematicsProperties(
      const std::vector<pid> &InitialState_,
      const std::vector<pid> &FinalState_,
      std::shared_ptr<PartList> ParticleList_,
      const ComPWA::FourMomentum &InitialStateP4_,
      const std::vector<unsigned int> &FinalStateEventPositionMapping_)
      : InitialState(InitialState_), FinalState(FinalState_),
        ParticleList(ParticleList_), InitialStateP4(InitialStateP4_),
        FinalStateEventPositionMapping(FinalStateEventPositionMapping_) {}
};

class Kinematics {
public:
  //! Constructor
  Kinematics(const KinematicsProperties &KinematicsProperties_);

  /// Delete copy constructor. For each Kinematics in the analysis only
  /// one instance should exist since Kinematics does the bookkeeping for which
  /// SubSystems variables needs to be calculated. That instance can then be
  /// passed as (smart) pointer. Note: Not sure if we also should delete the
  /// move constructor.
  Kinematics(const Kinematics &that) = delete;

  virtual ~Kinematics() {}

  /// Convert Event to DataPoint
  virtual void convert(const ComPWA::Event &ev, DataPoint &point) const = 0;

  /// Check if DataPoint is within phase space boundaries
  virtual bool isWithinPhsp(const DataPoint &point) const = 0;

  virtual double phspVolume() const;

  virtual void setPhspVolume(double phsp);

  virtual std::size_t numVariables() const { return VariableNames.size(); }

  /*virtual std::vector<pid> getFinalState() const { return KinematicsProperties.FinalState; }

  virtual std::vector<pid> getInitialState() const { return KinematicsProperties.InitialState; }

  virtual ComPWA::FourMomentum getInitialStateFourMomentum() const {
    return KinematicsProperties.InitialStateP4;
  }*/
  const ComPWA::KinematicsProperties& getKinematicsProperties() const {
	  return KinematicsProperties;
  }

  virtual unsigned int getDataID(const ComPWA::SubSystem &sys) const = 0;

  virtual std::vector<std::string> variableNames() const {
    return VariableNames;
  }

protected:
  virtual unsigned int
  convertFinalStateIDToPositionIndex(unsigned int fs_id) const;
  virtual std::vector<unsigned int> convertFinalStateIDToPositionIndex(
      const std::vector<unsigned int> &fs_ids) const;

  ComPWA::KinematicsProperties KinematicsProperties;

  /// Names of variables
  std::vector<std::string> VariableNames;

  virtual double calculatePhspVolume() const = 0;
  bool HasPhspVolume;
  double PhspVolume;
};

} // namespace ComPWA
#endif
