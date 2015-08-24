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

#ifndef PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_
#define PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_

#include "TopologyAmplitude.hpp"

#include "Core/Amplitude.hpp"
#include "Core/Efficiency.hpp"

namespace HelicityFormalism {

class CoherentAmplitude: public Amplitude {
  std::vector<TopologyAmplitude> topology_amplitudes_;
  //std::vector<HelicityDecayTree> decay_trees_;
  std::shared_ptr<Efficiency> efficiency_;

  ParameterList parameters_;
  std::vector<std::vector<IndexList> > data_point_index_lists_;

public:
  CoherentAmplitude(const std::vector<TopologyAmplitude>& amplitude_trees);
  virtual ~CoherentAmplitude();

  void registerTopologyAmplitudeParameters();

  void init(const Event& event);

  const double integral();
  const double integral(ParameterList& par);
  const double normalization();
  const double normalization(ParameterList& par);
  double getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen);
  double getMaxVal(std::shared_ptr<Generator> gen);

  const ParameterList& intensity(const dataPoint& point, ParameterList& par);
  const ParameterList& intensity(const dataPoint& point);
  const ParameterList& intensityNoEff(const dataPoint& point);
  const ParameterList& intensity(std::vector<double> point, ParameterList& par);

  const bool fillStartParVec(ParameterList& outPar);
  void setParameterList(ParameterList& par);

  void printAmps();
  void printFractions();

  Amplitude* Clone();
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_COHERENTAMPLITUDE_HPP_ */
