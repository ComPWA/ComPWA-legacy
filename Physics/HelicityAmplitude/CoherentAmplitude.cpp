//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//--------------------------------------------------------------------------------

#include "CoherentAmplitude.hpp"
#include "HelicityKinematics.hpp"


namespace HelicityFormalism {

CoherentAmplitude::CoherentAmplitude(
    const std::vector<TopologyAmplitude>& amplitude_trees) :
    topology_amplitudes_(amplitude_trees) {

}

CoherentAmplitude::~CoherentAmplitude() {
}

void CoherentAmplitude::init(const Event& event) {
  // TODO: I have to cast the kinematics instance to a helicity kinematics.
  // This is pretty bad practice... I really don't get the point behind the
  // singleton kinematics class...
  // If you do not instanciate a different type of kinematics class beforehand
  // you are fine, but otherwise we are in deep shit.
  HelicityKinematics* kinematics =
      (HelicityKinematics*) HelicityKinematics::createInstance();

  // initialize the kinematics class first
  kinematics->init(event);
  // now get the index lists that tell the topology amplitudes which
  // data point variables to use in their evaluation
  data_point_index_lists_ =
      kinematics->getTopologyAmplitudeDataPointIndexLists();
}

const double CoherentAmplitude::integral() {

}
const double CoherentAmplitude::integral(ParameterList& par) {

}
const double CoherentAmplitude::normalization() {

}
const double CoherentAmplitude::normalization(ParameterList& par) {

}
double CoherentAmplitude::getMaxVal(ParameterList& par,
    std::shared_ptr<Generator> gen) {

}
double CoherentAmplitude::getMaxVal(std::shared_ptr<Generator> gen) {

}

const ParameterList& CoherentAmplitude::intensity(std::vector<double> point,
    ParameterList& par) {
  setParameterList(par);
  dataPoint dataP(point);
  return intensity(dataP);
}

const ParameterList& CoherentAmplitude::intensity(const dataPoint& point,
    ParameterList& par) {
  setParameterList(par);
  return intensity(point);
}

const ParameterList& CoherentAmplitude::intensityNoEff(const dataPoint& point) {
  std::complex<double> intensity = 0;
  HelicityKinematics::instance();
  if (Kinematics::instance()->isWithinPhsp(point)) {

    // at first we have to create a vector of the kinematic variables
    // that we pass to the topology amplitudes
    // here the ordering has to be identical to the one in the data point
    // since the index lists refer to variable positions in the data point
    std::vector<KinematicVariables> kinematic_variables =
        createKinematicVariablesFromDataPoint(point);
    for (unsigned int i = 0; i < topology_amplitudes_.size(); ++i) {
      // get the correct helicity angles from the dataPoint
      const std::vector<IndexList> &topology_data_index_lists =
          data_point_index_lists_[i];
      for (unsigned int final_state_combination_index = 0;
          final_state_combination_index < topology_data_index_lists.size();
          ++final_state_combination_index) {

        intensity += topology_amplitudes_[i].evaluate(kinematic_variables,
            topology_data_index_lists[final_state_combination_index]);
      }
    }
    intensity = pow(std::abs(intensity), 2.0);
  }

  if (intensity != intensity) {
    BOOST_LOG_TRIVIAL(error)<<"Intensity is not a number!!";
    intensity = 0;
  }
  result.SetParameterValue(0, intensity);
  return result;
}

std::vector<KinematicVariables> CoherentAmplitude::createKinematicVariablesFromDataPoint(
    const dataPoint& point) const {
  std::vector<KinematicVariables> kinematic_variables;

  return kinematic_variables;
}

const ParameterList& CoherentAmplitude::intensity(const dataPoint& point) {
  intensityNoEff(point);
  double eff = efficiency_->evaluate(point);
  result.SetParameterValue(0, result.GetDoubleParameter(0)->GetValue() * eff);
  return result;
}

const bool CoherentAmplitude::fillStartParVec(ParameterList& outPar) {
  outPar = ParameterList(params_);
  return true;
}

void CoherentAmplitude::setParameterList(ParameterList& par) {
  //parameters varied by Minimization algorithm
  if (par.GetNDouble() != params_.GetNDouble())
    throw std::runtime_error(
        "setParameterList(): size of parameter lists don't match");
  for (unsigned int i = 0; i < params_.GetNDouble(); i++) {
    std::shared_ptr<DoubleParameter> p = params_.GetDoubleParameter(i);
    if (!p->IsFixed()) {
      p->SetValue(par.GetDoubleParameter(i)->GetValue());
      p->SetError(par.GetDoubleParameter(i)->GetError());
    }
  }
  return;
}

void CoherentAmplitude::printAmps() {

}
void CoherentAmplitude::printFractions() {

}

Amplitude* CoherentAmplitude::Clone() {
}

} /* namespace HelicityFormalism */
