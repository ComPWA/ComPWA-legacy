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
namespace MinLogLH {

MinLogLH::MinLogLH(std::shared_ptr<Amplitude> amp,
                   std::shared_ptr<DataReader::Data> data,
                   std::shared_ptr<DataReader::Data> phspSample,
                   std::shared_ptr<DataReader::Data> accSample,
                   unsigned int startEvent, unsigned int nEvents)
    : _useFunctionTree(0), nEvts_(0), nPhsp_(0), nStartEvt_(startEvent),
      nUseEvt_(nEvents), _dataSample(data), _phspSample(phspSample),
      _phspAccSample(accSample), _penaltyLambda(0.0) {
  _ampVec.push_back(amp);
  Init();

  return;
}

MinLogLH::MinLogLH(std::vector<std::shared_ptr<Amplitude>> ampVec,
                   std::vector<double> fraction,
                   std::shared_ptr<DataReader::Data> data,
                   std::shared_ptr<DataReader::Data> phspSample,
                   std::shared_ptr<DataReader::Data> accSample,
                   unsigned int startEvent, unsigned int nEvents)
    : _useFunctionTree(0), nEvts_(0), nPhsp_(0), nStartEvt_(startEvent),
      nUseEvt_(nEvents), _ampVec(ampVec), _fraction(fraction),
      _dataSample(data), _phspSample(phspSample), _phspAccSample(accSample),
      _penaltyLambda(0.0) {
  Init();

  return;
}

void MinLogLH::Init() {
  double sumFraction = std::accumulate(_fraction.begin(), _fraction.end(), 0.0);
  if (sumFraction > 1.0)
    throw std::runtime_error("MinLogLH::init() | Fractions sum larger 1.0!");

  if (_fraction.size() == _ampVec.size() - 1)
    _fraction.push_back(1 - sumFraction);

  // check size
  if (_fraction.size() != _ampVec.size())
    throw std::runtime_error("MinLogLH::init() | List of fractions "
                             "(" +
                             std::to_string(_fraction.size()) +
                             ")"
                             " does not match with list of amplitudes"
                             "(" +
                             std::to_string(_ampVec.size()) + ")"
                                                              "!");

  nEvts_ = _dataSample->getNEvents();
  nPhsp_ = _phspSample->getNEvents();
  if (!nUseEvt_)
    nUseEvt_ = nEvts_ - nStartEvt_;
  if (!(nStartEvt_ + nUseEvt_ <= nEvts_))
    nUseEvt_ = nEvts_ - nStartEvt_;
  if (!(nStartEvt_ + nUseEvt_ <= nPhsp_))
    nUseEvt_ = nPhsp_ - nStartEvt_;

  // Get data as ParameterList
  _dataSampleList = _dataSample->getListOfData();
  _phspSampleList = _phspSample->getListOfData();
  if (_phspAccSample)
    _phspAccSampleList = _phspAccSample->getListOfData();
  else
    _phspAccSampleList = _phspSample->getListOfData();

  calcSumOfWeights();

  calls = 0; // member of ControlParameter
}

void MinLogLH::Reset() {
  _ampVec.clear();
  _fraction.clear();
  _dataSample = std::shared_ptr<DataReader::Data>();
  _phspSample = std::shared_ptr<DataReader::Data>();
  _phspAccSample = std::shared_ptr<DataReader::Data>();
  _phspAccSampleEff = 1.0;
}

std::shared_ptr<ComPWA::ControlParameter>
MinLogLH::createInstance(std::shared_ptr<Amplitude> amp,
                         std::shared_ptr<DataReader::Data> data,
                         std::shared_ptr<DataReader::Data> phspSample,
                         unsigned int startEvent, unsigned int nEvents) {
  if (!instance_) {
    std::shared_ptr<DataReader::Data> accSample_ =
        std::shared_ptr<DataReader::Data>();
    instance_ = std::shared_ptr<ComPWA::ControlParameter>(
        new MinLogLH(amp, data, phspSample,
                     std::shared_ptr<DataReader::Data>(), // empty sample
                     startEvent, nEvents));
    LOG(debug) << "MinLogLH::createInstance() | "
                  "Creating instance from amplitude and dataset!";
  }
  return instance_;
}

std::shared_ptr<ComPWA::ControlParameter>
MinLogLH::createInstance(std::shared_ptr<Amplitude> amp,
                         std::shared_ptr<DataReader::Data> data,
                         std::shared_ptr<DataReader::Data> phspSample,
                         std::shared_ptr<DataReader::Data> accSample,
                         unsigned int startEvent, unsigned int nEvents) {
  if (!instance_) {
    instance_ = std::shared_ptr<ControlParameter>(
        new MinLogLH(amp, data, phspSample, accSample, startEvent, nEvents));
  }
  return instance_;
}

std::shared_ptr<ComPWA::ControlParameter>
MinLogLH::createInstance(std::vector<std::shared_ptr<Amplitude>> ampVec,
                         std::vector<double> frac,
                         std::shared_ptr<DataReader::Data> data_,
                         std::shared_ptr<DataReader::Data> phspSample_,
                         std::shared_ptr<DataReader::Data> accSample_,
                         unsigned int startEvent, unsigned int nEvents) {
  if (!instance_) {
    instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(
        ampVec, frac, data_, phspSample_, accSample_, startEvent, nEvents));
  }
  return instance_;
}

void MinLogLH::setAmplitude(std::shared_ptr<Amplitude> amp,
                            std::shared_ptr<DataReader::Data> data,
                            std::shared_ptr<DataReader::Data> phspSample,
                            std::shared_ptr<DataReader::Data> accSample,
                            unsigned int startEvent, unsigned int nEvents,
                            bool useFuncTr) {
  Reset();

  _ampVec.push_back(amp);
  _dataSample = data;
  _phspSample = phspSample;
  _phspAccSample = accSample;

  nStartEvt_ = startEvent;
  nUseEvt_ = nEvents;

  Init();

  return;
}

void MinLogLH::setAmplitude(std::vector<std::shared_ptr<Amplitude>> ampVec,
                            std::vector<double> frac,
                            std::shared_ptr<DataReader::Data> data,
                            std::shared_ptr<DataReader::Data> phspSample,
                            std::shared_ptr<DataReader::Data> accSample,
                            unsigned int startEvent, unsigned int nEvents,
                            bool useFuncTr) {
  _ampVec.clear();

  _ampVec = ampVec;
  _dataSample = data;
  _phspSample = phspSample;
  _phspAccSample = accSample;

  nStartEvt_ = startEvent;
  nUseEvt_ = nEvents;

  Init();

  return;
}

void MinLogLH::calcSumOfWeights() {
  _sumOfWeights = 0;
  for (unsigned int evt = nStartEvt_; evt < nUseEvt_ + nStartEvt_; evt++) {
    Event ev(_dataSample->getEvent(evt));
    _sumOfWeights += ev.getWeight();
  }
  LOG(info) << "MinLogLH: for current data set: numEvents = " << nUseEvt_
            << " sumOfWeights=" << _sumOfWeights << " for current data set.";
  return;
}

void MinLogLH::setUseFunctionTree(bool t) {
  if (t == 0)
    _useFunctionTree = 0;
  else
    iniLHtree();
}

void MinLogLH::iniLHtree() {
  LOG(debug) << "MinLogLH::iniLHtree() constructing the LH tree";

  if (_useFunctionTree)
    return;
  auto it = _ampVec.begin();
  for (; it != _ampVec.end(); ++it)
    if (!(*it)->hasTree())
      throw std::runtime_error("MinLogLH::iniLHtree() amplitude has no tree");

  //----Strategies needed
  std::shared_ptr<Strategy> mmultDStrat(new MultAll(ParType::MDOUBLE));
  std::shared_ptr<Strategy> multiDoubleAddStrat(new AddAll(ParType::MDOUBLE));
  std::shared_ptr<Strategy> mlogStrat(new LogOf(ParType::MDOUBLE));
  std::shared_ptr<Strategy> multDStrat(new MultAll(ParType::DOUBLE));
  std::shared_ptr<Strategy> addStrat(new AddAll(ParType::DOUBLE));
  std::shared_ptr<Strategy> invStrat(new Inverse(ParType::DOUBLE));

  LOG(debug) << "MinLogLH::iniLHTree() construction LH tree";
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

  unsigned int weightId = Kinematics::instance()->GetNVars() + 1;
  std::shared_ptr<MultiDouble> weight =
      _dataSampleList.GetMultiDouble(weightId);

  //-log L = (-1)*N/(\sum_{ev} w_{ev}) \sum_{ev} ...
  _tree->createHead("LH", multDStrat);
  _tree->createLeaf("minusOne", -1, "LH");
  _tree->createLeaf("nEvents", sampleSize, "LH");
  _tree->createNode("invSumWeights", invStrat, "LH");
  _tree->createNode("sumEvents", addStrat, "LH");
  _tree->createNode("sumWeights", addStrat, "invSumWeights");
  _tree->createLeaf("weight", weight, "sumWeights");
  _tree->createNode("weightLog", mmultDStrat, "sumEvents", sampleSize,
                    false); // w_{ev} * log( I_{ev} )
  _tree->createLeaf("weight", weight, "weightLog");
  _tree->createNode("Log", mlogStrat, "weightLog", sampleSize, false);
  // I_{ev} = x_{ev} + (1-f_{bkg})
  _tree->createNode("Add", multiDoubleAddStrat, "Log", sampleSize, false);

  for (int i = 0; i < _ampVec.size(); i++) {
    std::string name = _ampVec.at(i)->GetName();
    _tree->createNode("Intens_" + name, mmultDStrat, "Add", sampleSize, false);
    _tree->createLeaf("frac_" + name, _fraction.at(i), "Intens_" + name);
    // Expect that intensity is calculated and normalised within AmpAbsDyn
    _tree->insertTree(_ampVec.at(i)->GetTree(
                          _dataSampleList, _phspAccSampleList, _phspSampleList),
                      "Intens_" + name);
  }

  _tree->recalculate();
  if (!_tree->sanityCheck()) {
    throw std::runtime_error("MinLogLH::iniLHtree() | Tree has structural "
                             "problems. Sanity check not passed!");
  }
  LOG(debug) << "MinLogLH::iniLHtree() | "
                "Construction of LH tree finished!";
  _useFunctionTree = 1;
  return;
}

double MinLogLH::controlParameter(ParameterList &minPar) {
  double lh = 0;
  if (!_useFunctionTree) {
    // Calculate normalization
    double vol = Kinematics::instance()->GetPhspVolume();
    std::vector<double> normVec(_ampVec.size(), 0.0);

    std::shared_ptr<DataReader::Data> sam;
    if (_phspAccSample)
      sam = _phspAccSample;
    else
      sam = _phspSample;

    int size = sam->getNEvents();
    for (unsigned int phsp = 0; phsp < size; phsp++) { // loop over phspSample
      Event ev(sam->getEvent(phsp));
      dataPoint point;
      try {
        point = dataPoint(ev);
      } catch (BeyondPhsp &ex) { // event outside phase, remove
        continue;
      }
      for (unsigned int i = 0; i < _ampVec.size(); ++i) {
        ParameterList intensList = _ampVec.at(i)->intensity(point);
        double intens = intensList.GetDoubleParameter(0)->GetValue();
        if (intens > 0)
          normVec.at(i) += intens;
      }
    }
    for (auto it = normVec.begin(); it != normVec.end(); ++it) {
      if (!(*it)) {
        std::cout << "Amps: " << _ampVec.size() << " Norms: " << normVec.size()
                  << " Norm " << normVec.at(0) << std::endl;
        throw std::runtime_error(
            "MinLogLH::controlParameter() | "
            "Normalization can not be calculated for at least one "
            "resonance ");
      }
      (*it) *= vol / (double)size;
    }

    /* Approximate error of integration, see (Numerical Recipes Vol3,
     * p398, Eq. 7.7.1)
     * double normError = sqrt(
     * 			( vol*vol*normSq/sam_size - norm*norm) / sam_size
     * 		);
     */

    // Use internal amplitude integration - no unbinned eff correction
    // norm = amp->normalization();

    // Calculate \Sum_{ev} log()
    double sumLog = 0;
    // loop over data sample
    for (unsigned int evt = nStartEvt_; evt < nUseEvt_ + nStartEvt_; evt++) {
      Event ev(_dataSample->getEvent(evt));
      dataPoint point(ev);
      double val = 0;
      for (unsigned int i = 0; i < _ampVec.size(); ++i) {
        ParameterList intensList = _ampVec.at(i)->intensityNoEff(point);
        double tmp = intensList.GetDoubleParameter(0)->GetValue();
        val += tmp * _fraction.at(i) / normVec.at(i);
      }
      if (val > 0)
        sumLog += std::log(val) * ev.getWeight();
    }
    lh = (-1) * ((double)nUseEvt_) / _sumOfWeights * sumLog;
  } else {
    _tree->recalculate();
    std::shared_ptr<DoubleParameter> logLH =
        std::dynamic_pointer_cast<DoubleParameter>(_tree->head()->getValue());
    lh = logLH->GetValue();
  }
  lh += calcPenalty();
  calls++;
  return lh; // return -logLH
}

void MinLogLH::setPenaltyScale(double sc, int ampID) {
  if (sc < 0) {
    LOG(info) << "MinLogLH::setPenaltyScale | "
                 "Penalty scale cannot be negative!";
    return;
  }
  if (ampID < 0 || ampID >= _ampVec.size())
    LOG(info) << "MinLogLH::setPenaltyScale | "
                 "Amplitude ID not valid!";

  _penaltyAmpID = ampID;
  _penaltyLambda = sc;

  LOG(info) << "MinLogLH::setPenaltyScale | "
               "Setting scale of penalty term to "
            << sc << " for amplitude " << ampID << "!";
}

double MinLogLH::calcPenalty() {
  if (_penaltyLambda <= 0)
    return 0; // penalty term disabled
  double magSum = 0;
  auto amp = _ampVec.at(_penaltyAmpID);
  auto it = amp->GetResonanceItrFirst();
  std::vector<resonanceItr> resoList;
  for (; it != amp->GetResonanceItrLast(); ++it) {
    if ((*it)->GetName().find("_CP") != std::string::npos)
      continue;
    double v = std::fabs((*it)->GetMagnitude());
    magSum += v;
    resoList.push_back(it);
  }
  return (_penaltyLambda * magSum / std::sqrt(amp->GetIntegral(resoList)));
  //	return (_penaltyLambda*magSum);
}

// double MinLogLH::calcPenalty()
//{
//	if(_penaltyLambda<=0) return 0; //penalty term disabled
//	double ffSum = 0;
//	auto amp = _ampVec.at(_penaltyAmpID);
//	ParameterList ffList;
//	amp->GetFitFractions(ffList);
//
//	for(unsigned int i=0;i<ffList.GetNDouble(); ++i){
//		double norm = amp->GetIntegral(resoList);
//
//		ffSum += ffList.GetDoubleParameter(i)->GetValue();
//	}
//	return (_penaltyLambda*ffSum);
//}

} /* namespace MinLogLH */
} /* namespace Estimator */
} /* namespace ComPWA */
