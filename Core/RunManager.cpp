//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding functionality to generate set of
//events
//-------------------------------------------------------------------------------
#include <memory>
#include <ctime>
#include <numeric>

#include <time.h>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/progress.hpp>

#include "DataReader/Data.hpp"
#include "Estimator/Estimator.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Event.hpp"
#include "Core/Generator.hpp"
#include "Core/ProgressBar.hpp"

#include "Core/RunManager.hpp"

namespace ComPWA {

using namespace boost::log;

RunManager::RunManager() {}

RunManager::RunManager(std::shared_ptr<DataReader::Data> data,
                       std::shared_ptr<AmpIntensity> intens,
                       std::shared_ptr<Optimizer::Optimizer> optimizer)
    : sampleData_(data), opti_(optimizer), intens_(intens){}

RunManager::RunManager(unsigned int size, std::shared_ptr<AmpIntensity> intens,
                       std::shared_ptr<Generator> gen)
    :  gen_(gen), intens_(intens) {}

RunManager::~RunManager() {
  if( gen_ )
    LOG(debug) << "~RunManager: Last seed: " << gen_->GetSeed();
}

std::shared_ptr<FitResult> RunManager::Fit(ParameterList &inPar) {
  LOG(info) << "RunManager::startFit() | Starting minimization.";

  // MINIMIZATION
  std::shared_ptr<FitResult> result = opti_->exec(inPar);

  LOG(info) << "RunManager::startFit() | Minimization finished!"
               " Result = "
            << result->GetResult() << ".";

  return result;
}

void RunManager::SetPhspSample(std::shared_ptr<Data> phsp,
                               std::shared_ptr<DataReader::Data> truePhsp) {
  if (truePhsp && truePhsp->getNEvents() != phsp->getNEvents())
    throw std::runtime_error(
        "RunManager::setPhspSample() | "
        "Reconstructed sample and true sample have not the same size!");
  samplePhsp_ = phsp;
  sampleTruePhsp_ = truePhsp;
}

void RunManager::SetTruePhspSample(std::shared_ptr<DataReader::Data> truePhsp) {
  if (truePhsp && samplePhsp_ &&
      truePhsp->getNEvents() != samplePhsp_->getNEvents())
    throw std::runtime_error(
        "RunManager::setPhspSample() | "
        "Reconstructed sample and true sample have not the same size!");

  sampleTruePhsp_ = truePhsp;
}

bool RunManager::Generate(int number) {
  LOG(info) << "RunManager::generate() | "
               "Generating "
            << number << " signal events!";

  return gen(number, gen_, intens_, sampleData_, samplePhsp_, sampleTruePhsp_);
}

bool RunManager::gen(int number, std::shared_ptr<Generator> gen,
                     std::shared_ptr<AmpIntensity> amp,
                     std::shared_ptr<DataReader::Data> data,
                     std::shared_ptr<DataReader::Data> phsp,
                     std::shared_ptr<DataReader::Data> phspTrue) {

  if (number == 0)
    return 0;

  // Doing some checks
  if (number < 0 && !phsp)
    throw std::runtime_error("RunManager: gen() negative number of events: " +
                             std::to_string((long double)number) +
                             ". And no phsp sample given!");
  if (!amp)
    throw std::runtime_error("RunManager::gen() | Amplitude not valid");
  if (!gen)
    throw std::runtime_error("RunManager::gen() | Generator not valid");
  if (!data)
    throw std::runtime_error("RunManager::gen() | Sample not valid");
  if (data->getNEvents() > 0)
    throw std::runtime_error("RunManager::gen() | Sample not empty!");
  if (phspTrue && !phsp)
    throw std::runtime_error("RunManager::gen() | We have a sample of true"
                             " phsp events, but no phsp sample!");
  if (phspTrue && phspTrue->getNEvents() != phsp->getNEvents())
    throw std::runtime_error(
        "RunManager::gen() | We have a sample of true "
        "phsp events, but the sample size doesn't match that one of "
        "the phsp sample!");

  double maxSampleWeight = 1.0;
  if (phsp)
    maxSampleWeight = phsp->getMaxWeight();
  if (phspTrue && phspTrue->getMaxWeight() > maxSampleWeight)
    maxSampleWeight = phspTrue->getMaxWeight();

  /* Maximum value for random number generation. We introduce an arbitrary
   * factor of 1.2 to make sure that the maximum value is never reached.
   */
  double generationMaxValue = 5 * amp->GetMaximum(gen) * maxSampleWeight;

  unsigned int initialSeed = gen->GetSeed();
  unsigned int totalCalls = 0;
  //	unsigned int startTime = clock();
  double AMPpdf;

  unsigned int limit;
  unsigned int acceptedEvents = 0;
  if (phsp) {
    limit = phsp->getNEvents();
  } else
    limit = 100000000; // set large limit, should never be reached

  Event evt;     // event that we fill into generated sample
  Event evtTrue; // event that is used to evalutate amplitude
  progressBar bar(number);
  if (number <= 0)
    bar = progressBar(limit);
  for (unsigned int i = 0; i < limit; i++) {
    if (phsp && phspTrue) { // phsp and true sample is set
      evtTrue = phspTrue->getEvent(i);
      evt = phsp->getEvent(i);
    } else if (phsp) { // phsp sample is set
      evt = phsp->getEvent(i);
      evtTrue = evt;
    } else { // otherwise generate event
      gen->Generate(evt);
      evtTrue = evt;
    }
    if (number <= 0)
      bar.nextEvent();

    // use reconstructed position for weights
    double weight = evt.getWeight();

    // use true position for amplitude value
    dataPoint point;
    try {
      point = dataPoint(evtTrue);
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }

    totalCalls++;
    double ampRnd = gen->GetUniform(0, generationMaxValue);
    AMPpdf = amp->Intensity(point);

    if (generationMaxValue < (AMPpdf * weight)){
      LOG(error) << "RunManager::gen() | Error in HitMiss "
          "procedure: Maximum value of random number generation "
          "smaller then amplitude maximum! We raise the maximum value and restart generation!";
      i=0;
      bar = progressBar(number);
      gen->SetSeed(initialSeed);
      generationMaxValue = 2 * (AMPpdf * weight);
      data->Clear();
      continue;
//      throw std::runtime_error(
//          "RunManager::gen() | Error in HitMiss "
//          "procedure: Maximum value of random number generation "
//          "smaller then amplitude maximum!");
    }

    if (ampRnd > (weight * AMPpdf))
      continue;

    // Fill event to sample
    /* reset weights: the weights are taken into account by hit and miss. The
     * resulting
     * sample is therefore unweighted */
    evt.setWeight(1.);     // reset weight
    evt.setEfficiency(1.); // reset weight
    data->pushEvent(evt);  // unfortunatly not thread safe

    // some statistics
    acceptedEvents++;
    if (number > 0)
      bar.nextEvent();
    // break if we have a sufficienct number of events
    if (acceptedEvents >= number)
      i = limit;
  }
  if (data->getNEvents() < number)
    LOG(error) << "RunManager::gen() | Not able to generate "
                  " "
               << number << " events. Phsp sample too small. Current size "
                            "of sample is now "
               << data->getNEvents();

  LOG(info) << "Efficiency of toy MC generation: "
            << (double)data->getNEvents() / totalCalls;

  return true;
}

bool RunManager::GeneratePhsp(int number) {
  if (number == 0)
    return 0;
  if (!samplePhsp_)
    throw std::runtime_error("RunManager: generatePhsp() | "
                             "No phase-space sample set");
  if (samplePhsp_->getNEvents() > 0)
    throw std::runtime_error("RunManager: generatePhsp() | "
                             "Dataset not empty! abort!");

  LOG(info) << "Generating phase-space MC: [" << number << " events] ";

  progressBar bar(number);
  for (unsigned int i = 0; i < number; i++) {
    if (i > 0)
      i--;
    Event tmp;
    gen_->Generate(tmp);
    double ampRnd = gen_->GetUniform(0,1);
    if (ampRnd > tmp.getWeight())
      continue;

    /* Reset weights: weights are taken into account by hit&miss. The
     * resulting sample is therefore unweighted
     */
    tmp.setWeight(1.);

    tmp.setEfficiency(1.);
    i++;
    samplePhsp_->pushEvent(tmp); // unfortunatly not thread safe
    bar.nextEvent();
  }
  return true;
}

} /* namespace ComPWA */
