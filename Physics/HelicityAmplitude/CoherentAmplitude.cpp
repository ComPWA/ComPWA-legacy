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
    topology_amplitudes_(amplitude_trees), wasMaxAmplitudeValueCalculated_(
        false), maxAmplitudeValue_(0.0) {

  registerTopologyAmplitudeParameters();
}

CoherentAmplitude::~CoherentAmplitude() {
}

void CoherentAmplitude::registerTopologyAmplitudeParameters() {
  for (unsigned int topology_amplitude_index = 0;
      topology_amplitude_index < topology_amplitudes_.size();
      ++topology_amplitude_index) {

    const std::vector<SequentialTwoBodyDecayAmplitude>& sequential_decay_list =
        topology_amplitudes_[topology_amplitude_index].getSequentialDecayList();

    std::vector<SequentialTwoBodyDecayAmplitude>::const_iterator sequential_decay_iter;

    for (sequential_decay_iter = sequential_decay_list.begin();
        sequential_decay_iter != sequential_decay_list.end();
        ++sequential_decay_iter) {

      parameters_.AddParameter(sequential_decay_iter->strength_);
      parameters_.AddParameter(sequential_decay_iter->phase_);

      std::vector<FullTwoBodyDecayAmplitude>::const_iterator decay_node_iter;

      for (decay_node_iter =
          sequential_decay_iter->full_decay_amplitude_chain_list_.begin();
          decay_node_iter
              != sequential_decay_iter->full_decay_amplitude_chain_list_.end();
          ++decay_node_iter) {
        parameters_.Append(decay_node_iter->second->getParameterList());
      }
    }
  }
}

void CoherentAmplitude::init() {
// TODO: I have to cast the kinematics instance to a helicity kinematics.
// This is pretty bad practice... I really don't get the point behind the
// singleton kinematics class...
// If you do not instantiate a different type of kinematics class beforehand
// you are fine, but otherwise we are in deep shit.
  HelicityKinematics* kinematics =
      (HelicityKinematics*) HelicityKinematics::createInstance();

// initialize the kinematics class first
//  kinematics->init(event);
// now get the index lists that tell the topology amplitudes which
// data point variables to use in their evaluation
  data_point_index_lists_ =
      kinematics->getTopologyAmplitudeDataPointIndexLists();

  // and initialize parameter at position 0 in the result objects as result
  std::shared_ptr<DoubleParameter> result_value(
      new DoubleParameter("coherent_amp"));
  result.AddParameter(result_value);
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
  setParameterList(par);
  return getMaxVal(gen);
}
double CoherentAmplitude::getMaxVal(std::shared_ptr<Generator> gen) {
  if (!wasMaxAmplitudeValueCalculated_) {
    HelicityKinematics* kin =
        dynamic_cast<HelicityKinematics*>(Kinematics::instance());

    Event event;
    double maxVal = 0;
    for (unsigned int i = 0; i < 100000; i++) {
      //create event
      gen->generate(event);
      dataPoint point;
      kin->eventToDataPoint(event, point);

      if (!kin->isWithinPhsp(point)) {
        if (i > 0)
          i--;
        continue;
      }    //only integrate over phase space
      ParameterList res = intensity(point);
      double intens = *res.GetDoubleParameter(0);
      if (intens > maxVal) {
        maxVal = intens;
      }
    }

    maxAmplitudeValue_ = maxVal;
    wasMaxAmplitudeValueCalculated_ = true;

    BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity::calcMaxVal() calculated maximum of amplitude: "
    <<maxAmplitudeValue_;
  }

  return maxAmplitudeValue_;
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
  std::complex<double> coherent_amplitude = 0;
  double intensity(0.0);

  if (Kinematics::instance()->isWithinPhsp(point)) {
    for (unsigned int i = 0; i < topology_amplitudes_.size(); ++i) {
      // at first get the appropriate data evaluation index lists for this
      // topology which we pass to the topology amplitudes later
      const std::vector<IndexList> &topology_data_index_lists =
          data_point_index_lists_[i];
      for (unsigned int final_state_combination_index = 0;
          final_state_combination_index < topology_data_index_lists.size();
          ++final_state_combination_index) {
        coherent_amplitude += topology_amplitudes_[i].evaluate(point,
            topology_data_index_lists[final_state_combination_index]);
      }
    }
    intensity = pow(std::abs(coherent_amplitude), 2.0);
  }

  if (intensity != intensity) {
    BOOST_LOG_TRIVIAL(error)<<"Intensity is not a number!!";
    intensity = 0;
  }
  result.SetParameterValue(0, intensity);
  return result;
}

const ParameterList& CoherentAmplitude::intensity(const dataPoint& point) {
  intensityNoEff(point);
  if (0 != efficiency_.get()) {
    double eff = efficiency_->evaluate(point);
    result.SetParameterValue(0, result.GetDoubleParameter(0)->GetValue() * eff);
  }
  return result;
}

const bool CoherentAmplitude::fillStartParVec(ParameterList& outPar) {
  outPar = ParameterList(parameters_);
  return true;
}

void CoherentAmplitude::setParameterList(ParameterList& par) {
//parameters varied by Minimization algorithm
  if (par.GetNDouble() != parameters_.GetNDouble())
    throw std::runtime_error(
        "setParameterList(): size of parameter lists don't match");
  for (unsigned int i = 0; i < parameters_.GetNDouble(); i++) {
    std::shared_ptr<DoubleParameter> p = parameters_.GetDoubleParameter(i);
    if (!p->IsFixed()) {
      p->SetValue(par.GetDoubleParameter(i)->GetValue());
      p->SetError(par.GetDoubleParameter(i)->GetError());
    }
  }
  return;
}

bool CoherentAmplitude::copyParameterList(ParameterList& par) {
  par = parameters_;
  return true;
}

double CoherentAmplitude::getIntValue(std::string var1, double min1,
    double max1, std::string var2, double min2, double max2) {
  return 0.0;
}

void CoherentAmplitude::printAmps() {

}
void CoherentAmplitude::printFractions() {

}

Amplitude * CoherentAmplitude::Clone() {
}

} /* namespace HelicityFormalism */
