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

#ifndef PHYSICS_HELICITYAMPLITUDE_HELICITYKINEMATICS_HPP_
#define PHYSICS_HELICITYAMPLITUDE_HELICITYKINEMATICS_HPP_

#include "Core/Kinematics.hpp"

//#include "Physics/HelicityAmplitude/ParticleStateDefinitions.hpp"
//#include "Physics/HelicityAmplitude/FinalStateParticleCombinatorics.hpp"

#include "Physics/qft++/Vector4.h"

class Particle;

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class HelicityKinematics: public Kinematics {
  
public:
  
  static Kinematics* createInstance(std::vector<int> initialState, std::vector<int> finalStateIds) {
    if ( _inst ) {
      throw std::runtime_error("HelicityKinematics::createInstance() | Instance already exists!");
    }
    _inst = new HelicityKinematics(initialState, finalStateIds);
    return _inst;
  }
  
  //! converts Event to dataPoint
  void EventToDataPoint(const Event& event, dataPoint& point) const;
  
  // delete methods to ensure that there will only be one instance
  HelicityKinematics() = delete;
  
  HelicityKinematics(const HelicityKinematics&) = delete;
  
  void operator=(const HelicityKinematics&) = delete;

  bool IsWithinPhsp (const dataPoint& point) const;

  virtual bool IsWithinBoxPhsp(int idA, int idB, double varA,
                               double varB) const { return true; }
  
  //! get mass of particles
  virtual double GetMass(unsigned int num) const { return 0; }
  
  //! get mass of particles
  virtual double GetMass(std::string name) const { return 0; }
  
  //! Get name of mother particle
  virtual std::string GetMotherName() const { return _nameMother; };
  
  //! Get mass of mother particle
  virtual double GetMotherMass() const { return _M; }
  
  //! Get number of particles
  virtual unsigned int GetNumberOfParticles() const { return _finalState.size(); }
  
  virtual std::vector<int> GetDataID( int id1, int id2 ){
    return std::vector<int>(1,1);
  }
  
//  FinalStateParticleCombinatorics fsp_combinatorics_;
//  std::map<std::vector<IDInfo>, std::vector<IndexList> > unique_id_lists_;
//
//  // this has the decay information of the topology amplitudes
//  std::vector<TwoBodyDecayTopology> decay_topologies_;

  // ---------------------------------------------------------------------
  // the following containers are only filled once, during the
  // initialization step and are relevant for the data point calculations
  // ---------------------------------------------------------------------
  // one entry for each occuring cms combination
//  std::vector<IndexList> unique_occurring_decay_event_fs_index_lists_;
  // one entry for each occuring decay into two definite cms products
//  std::vector<TwoBodyDecayIndices> unique_occurring_decay_product_index_list_links_;

  // outer vector: one entry for each decay topology
  // inner vector: one entry for each final state particle grouping combination
  // (can be more than one for indistinguishable fs particles)
  // IndexList: list of indices to the value of the data point used in the
  // evaluation of the topology amplitudes
//  std::vector<std::vector<IndexList> > topology_amplitude_data_point_index_lists_;

  // kinematic variables index locations for one individual set
//  std::map<std::string, unsigned int> kinematic_variable_to_index_mapping;
  // ---------------------------------------------------------------------


//  void buildDataPointIndexListForTopology(unsigned int index,
//      const TwoBodyDecayTopology& topology,
//      const IndexMapping& fs_particle_mapping,
//      IndexList& data_point_index_list);
//
//  unsigned int convertAndStoreParticleIndexList(const IndexList& particle_list,
//      const IndexMapping& fs_particle_mapping);
//
//  IndexList convertParticleIDToEventIndexList(const IndexList& particle_id_list,
//      const IndexMapping& fs_particle_mapping) const;
//
  std::vector<Vector4<double> > createRequired4Vectors(
      const Event& event) const;

  void addParticleToCMS4Vector(const Particle& event_particle,
      Vector4<double>& cms_4vector) const;

  void fillPointWithBoostedKinematicVariables(dataPoint& point,
      const std::vector<Vector4<double> >& unique_occurring_cms_4vectors,
      const TwoBodyDecayIndices& two_body_state_indices,
      unsigned int &data_point_fill_position) const;

  double calculateDecayAngle(const Vector4<double> &p, const Vector4<double> &q,
      const Vector4<double> &d) const;

protected:
  HelicityKinematics(std::vector<int> initialState, std::vector<int> finalState);
  
  virtual ~HelicityKinematics();
  
  double calculatePSArea();

  std::vector<std::string> _finalState;
  
  // Parameters of decaying mother particle (we assume that we have a decay)
  std::string _nameMother; //! name of mother particle
  int _idMother;   //! name of mother particle
  double _M;               //! mass of mother particle
  double _Msq;             //! mass of mother particle
  ComPWA::Spin _spinM;     //! spin of mother particle

//  void setDecayTopologies(
//      const std::vector<TwoBodyDecayTopology>& decay_topologies);
//
//  void init(const FinalStateParticleCombinatorics& fsp_combinatorics);

//  std::vector<std::vector<IndexList> > getTopologyAmplitudeDataPointIndexLists() const;


};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_HELICITYKINEMATICS_HPP_ */
