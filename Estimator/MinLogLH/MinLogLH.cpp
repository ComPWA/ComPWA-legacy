//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//	   Peter Weidenkaff - Weights and background fractions
//-------------------------------------------------------------------------------
#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <exception>
#include <numeric>

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Kinematics.hpp"
#include "Core/FitResult.hpp"

namespace ComPWA {
namespace Estimator {

MinLogLH::MinLogLH(Kinematics* kin,
                   std::shared_ptr<AmpIntensity> intens,
                   std::shared_ptr<DataReader::Data> data,
                   std::shared_ptr<DataReader::Data> phspSample,
                   std::shared_ptr<DataReader::Data> accSample,
                   unsigned int startEvent, unsigned int nEvents)
    : kin_(kin), _intens(intens), nEvts_(0), nPhsp_(0), nStartEvt_(startEvent),
      nUseEvt_(nEvents), _dataSample(data), _phspSample(phspSample),
      _phspAccSample(accSample) {

  Init();

  return;
}

void MinLogLH::Init() {

  nPhsp_ = _phspSample->GetNEvents();
  if (!nUseEvt_)
    nUseEvt_ = nEvts_ - nStartEvt_;
  if (!(nStartEvt_ + nUseEvt_ <= nEvts_))
    nUseEvt_ = nEvts_ - nStartEvt_;
  if (!(nStartEvt_ + nUseEvt_ <= nPhsp_))
    nUseEvt_ = nPhsp_ - nStartEvt_;

  // Get data as ParameterList
  _dataSampleList = _dataSample->GetListOfData();
  _phspSampleList = _phspSample->GetListOfData();
  if (_phspAccSample)
    _phspAccSampleList = _phspAccSample->GetListOfData();
  else
    _phspAccSampleList = _phspSample->GetListOfData();

  CalcSumOfWeights();

  LOG(info) << "MinLogLH::Init() |  Size of data sample = " << nUseEvt_
            << " ( Sum of weights = " << _sumOfWeights << " ).";

  calls = 0; // member of ControlParameter
}

void MinLogLH::Reset() {
  _intens = std::shared_ptr<AmpIntensity>();
  _dataSample = std::shared_ptr<DataReader::Data>();
  _phspSample = std::shared_ptr<DataReader::Data>();
  _phspAccSample = std::shared_ptr<DataReader::Data>();
  _phspAccSampleEff = 1.0;
}

std::shared_ptr<ComPWA::ControlParameter>
MinLogLH::CreateInstance(Kinematics* kin,
                         std::shared_ptr<AmpIntensity> intens,
                         std::shared_ptr<DataReader::Data> data,
                         std::shared_ptr<DataReader::Data> phspSample,
                         unsigned int startEvent, unsigned int nEvents) {

  if (!instance_) {
    std::shared_ptr<DataReader::Data> accSample_ =
        std::shared_ptr<DataReader::Data>();
    instance_ = std::shared_ptr<ComPWA::ControlParameter>(
        new MinLogLH(kin, intens, data, phspSample,
                     std::shared_ptr<DataReader::Data>(), // empty sample
                     startEvent, nEvents));
    LOG(debug) << "MinLogLH::createInstance() | "
                  "Creating instance from amplitude and dataset!";
  }
  return instance_;
}

std::shared_ptr<ComPWA::ControlParameter>
MinLogLH::CreateInstance(Kinematics* kin,
                         std::shared_ptr<AmpIntensity> intens,
                         std::shared_ptr<DataReader::Data> data,
                         std::shared_ptr<DataReader::Data> phspSample,
                         std::shared_ptr<DataReader::Data> accSample,
                         unsigned int startEvent, unsigned int nEvents) {
  if (!instance_) {
    instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(
        kin, intens, data, phspSample, accSample, startEvent, nEvents));
  }
  return instance_;
}

void MinLogLH::CalcSumOfWeights() {
  _sumOfWeights = 0;
  for (unsigned int evt = nStartEvt_; evt < nUseEvt_ + nStartEvt_; evt++) {
    Event ev(_dataSample->GetEvent(evt));
    _sumOfWeights += ev.GetWeight();
  }
  return;
}

void MinLogLH::IniLHtree() {
  LOG(debug) << "MinLogLH::IniLHtree() | Constructing FunctionTree!";
  if (!_intens->HasTree())
    throw std::runtime_error("MinLogLH::IniLHtree() |  AmpIntensity does not "
                             "provide a FunctionTree!");

  /* CONSTRUCTION OF THE LIKELIHOOD:
   * We denote the coherent sum over all resonances with T:
   * 		T := \sum_{i,j} c_i c_j^*A_iA_j^*
   * The negative log LH is given by:
   * 		-log L = - N/(\sum_{ev} w_{ev}) \sum_{ev} w_{ev} \log{f_{bkg}
   * \frac{|T|^2}{\int_{DP} |T|^2} + (1-f_{bkg})}
   * The sum over all weights is necessary to normalize the weights to one.
   * Otherwise the error
   * estimate is incorrect. The LH normalization is norm_{LH} = \int_{DP} |T|^2.
   * This formulation includes event weights as well as a flat background
   * description. f_{bkg} is
   * the fraction of background in the sample. Using both is of course
   * non-sense. Set weights to
   * one OR f_{bkg} to zero.
   */
  _tree = std::shared_ptr<FunctionTree>(new FunctionTree());
  int sampleSize = _dataSampleList.GetMultiDouble(0)->GetNValues();

  //-log L = (-1)*N/(\sum_{ev} w_{ev}) \sum_{ev} ...
  _tree->createHead("LH",
                    std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)));
  _tree->createLeaf("minusOne", -1, "LH");
  _tree->createLeaf("nEvents", sampleSize, "LH");
  _tree->createNode("invSumWeights",
                    std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)),
                    "LH");
  _tree->createNode("sumEvents",
                    std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                    "LH");
  _tree->createLeaf("SumOfWeights", _sumOfWeights, "invSumWeights");
  _tree->createNode("weightLog",
                    std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)),
                    "sumEvents", sampleSize,
                    false); // w_{ev} * log( I_{ev} )
//  _tree->createLeaf("weight", weight, "weightLog");
  _tree->createNode("Log",
                    std::shared_ptr<Strategy>(new LogOf(ParType::MDOUBLE)),
                    "weightLog", sampleSize, false);
  _tree->insertTree(_intens->GetTree(_dataSampleList, _phspAccSampleList,
                                     _phspSampleList, kin_->GetNVars()),
                    "Log");

  _tree->recalculate();
  if (!_tree->sanityCheck()) {
    throw std::runtime_error("MinLogLH::IniLHtree() | Tree has structural "
                             "problems. Sanity check not passed!");
  }
  LOG(debug) << "MinLogLH::IniLHtree() | "
                "Construction of LH tree finished!";
  return;
}

double MinLogLH::controlParameter(ParameterList &minPar) {
  double lh = 0;
  if (!_tree) {
    // Calculate \Sum_{ev} log()
    double sumLog = 0;
    // loop over data sample
    for (unsigned int evt = nStartEvt_; evt < nUseEvt_ + nStartEvt_; evt++) {
      Event ev(_dataSample->GetEvent(evt));
      dataPoint point(ev);
      double val = _intens->Intensity(point);
      sumLog += std::log(val) * ev.GetWeight();
    }
    lh = (-1) * ((double)nUseEvt_) / _sumOfWeights * sumLog;
  } else {
    _tree->recalculate();
    std::shared_ptr<DoubleParameter> logLH =
        std::dynamic_pointer_cast<DoubleParameter>(_tree->head()->getValue());
    lh = logLH->GetValue();
  }
  //  lh += calcPenalty();
  calls++;
  return lh; // return -logLH
}

} /* namespace Estimator */
} /* namespace ComPWA */
