// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <memory>
#include <ctime>
#include <numeric>

#include <time.h>

#include <boost/progress.hpp>

#include "DataReader/Data.hpp"
#include "Core/Estimator.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Event.hpp"
#include "Core/Generator.hpp"
#include "Core/ProgressBar.hpp"
#include "Tools/Integration.hpp"

#include "Tools/RunManager.hpp"

namespace ComPWA {

using namespace boost::log;

RunManager::RunManager() {}

RunManager::RunManager(std::shared_ptr<DataReader::Data> data,
                       std::shared_ptr<AmpIntensity> intens,
                       std::shared_ptr<Optimizer::Optimizer> optimizer)
    : sampleData_(data), opti_(optimizer), intens_(intens) {}

RunManager::RunManager(unsigned int size, std::shared_ptr<AmpIntensity> intens,
                       std::shared_ptr<Generator> gen)
    : gen_(gen), intens_(intens) {}

RunManager::~RunManager() {
  if (gen_)
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
  if (truePhsp && truePhsp->GetNEvents() != phsp->GetNEvents())
    throw std::runtime_error(
        "RunManager::setPhspSample() | "
        "Reconstructed sample and true sample have not the same size!");
  samplePhsp_ = phsp;
  sampleTruePhsp_ = truePhsp;
}

void RunManager::SetTruePhspSample(std::shared_ptr<DataReader::Data> truePhsp) {
  if (truePhsp && samplePhsp_ &&
      truePhsp->GetNEvents() != samplePhsp_->GetNEvents())
    throw std::runtime_error(
        "RunManager::setPhspSample() | "
        "Reconstructed sample and true sample have not the same size!");

  sampleTruePhsp_ = truePhsp;
}

bool RunManager::Generate(std::shared_ptr<Kinematics> kin, int number) {
  LOG(info) << "RunManager::generate() | "
               "Generating "
            << number << " signal events!";

  return gen(number, kin, gen_, intens_, sampleData_, samplePhsp_, sampleTruePhsp_);
}

bool RunManager::gen(int number, std::shared_ptr<Kinematics> kin, std::shared_ptr<Generator> gen,
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
  if (data->GetNEvents() > 0)
    throw std::runtime_error("RunManager::gen() | Sample not empty!");
  if (phspTrue && !phsp)
    throw std::runtime_error("RunManager::gen() | We have a sample of true"
                             " phsp events, but no phsp sample!");
  if (phspTrue && phspTrue->GetNEvents() != phsp->GetNEvents())
    throw std::runtime_error(
        "RunManager::gen() | We have a sample of true "
        "phsp events, but the sample size doesn't match that one of "
        "the phsp sample!");

  double maxSampleWeight = 1.0;
  if (phsp)
    maxSampleWeight = phsp->GetMaxWeight();
  if (phspTrue && phspTrue->GetMaxWeight() > maxSampleWeight)
    maxSampleWeight = phspTrue->GetMaxWeight();

  // Maximum value for random number generation. We introduce an arbitrary
  // factor of 5 to make sure that the maximum value is never reached.
  double generationMaxValue = 5* maxSampleWeight;

  unsigned int initialSeed = gen->GetSeed();
  unsigned int totalCalls = 0;

  unsigned int limit;
  if (phsp) {
    limit = phsp->GetNEvents();
    generationMaxValue *= Tools::Maximum(kin, amp, phsp) ;
  } else {
    limit = 100000000; // set large limit, should never be reached
  }

  LOG(trace) << "RunMananger::gen() | Using " << generationMaxValue
             << " as maximum value of the intensity.";
  
  Event evt;     // event that we fill into generated sample
  Event evtTrue; // event that is used to evalutate amplitude
  progressBar bar(number);
  if (number <= 0)
    bar = progressBar(limit);
  for (unsigned int i = 0; i < limit; i++) {
    if (phsp && phspTrue) { // phsp and true sample is set
      evtTrue = phspTrue->GetEvent(i);
      evt = phsp->GetEvent(i);
    } else if (phsp) { // phsp sample is set
      evt = phsp->GetEvent(i);
      evtTrue = evt;
    } else { // otherwise generate event
      gen->Generate(evt);
      evtTrue = evt;
    }
    if (number <= 0)
      bar.nextEvent();

    // use reconstructed position for weights
    double weight = evt.GetWeight();

    // use true position for amplitude value
    dataPoint point;
    try {
      kin->EventToDataPoint(evtTrue, point);
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }

    totalCalls++;
    double ampRnd = gen->GetUniform(0, generationMaxValue);
    double AMPpdf = amp->Intensity(point);

    // If maximum of intensity is reached we have to restart the procedure
    if (generationMaxValue < (AMPpdf * weight)) {
      LOG(trace) << "RunManager::gen() | Error in HitMiss "
                    "procedure: Maximum value of random number generation "
                    "smaller then amplitude maximum! We raise the maximum "
                    "value and restart generation!";
      i = 0;
      bar = progressBar(number);
      gen->SetSeed(initialSeed);
      generationMaxValue = 2 * (AMPpdf * weight);
      data->Clear();
      totalCalls = 0;
      continue;
    }

    if (ampRnd > (weight * AMPpdf))
      continue;

    // Fill event to sample
    // reset weights: the weights are taken into account by hit and miss. The
    // resulting sample is therefore unweighted
    evt.SetWeight(1.);     // reset weight
    evt.SetEfficiency(1.); // reset weight
    data->PushEvent(evt);

    if (number > 0)
      bar.nextEvent();
    // break if we have a sufficienct number of events
    if (data->GetNEvents() >= number)
      i = limit;
  }
  if (data->GetNEvents() < number)
    LOG(error) << "RunManager::gen() | Not able to generate "
               << number << " events. Phsp sample too small. Current size "
                            "of sample is now " << data->GetNEvents();

  LOG(info) << "Efficiency of toy MC generation: "
            << (double)data->GetNEvents() / totalCalls;

  return true;
}

bool RunManager::GeneratePhsp(int number) {
  if (number == 0)
    return 0;
  if (!samplePhsp_)
    throw std::runtime_error("RunManager: generatePhsp() | "
                             "No phase-space sample set");
  if (samplePhsp_->GetNEvents() > 0)
    throw std::runtime_error("RunManager: generatePhsp() | "
                             "Dataset not empty! abort!");

  LOG(info) << "Generating phase-space MC: [" << number << " events] ";

 boost::progress_display show_progress(number);
//  progressBar bar(number);
  for (unsigned int i = 0; i < number; i++) {
    if (i > 0)
      i--;
    Event tmp;
    gen_->Generate(tmp);
    double ampRnd = gen_->GetUniform(0, 1);
    if (ampRnd > tmp.GetWeight())
      continue;

    // Reset weights: weights are taken into account by hit&miss. The
    // resulting sample is therefore unweighted
    tmp.SetWeight(1.);

    tmp.SetEfficiency(1.);
    i++;
    samplePhsp_->PushEvent(tmp); // unfortunatly not thread safe
//    bar.nextEvent();
	  ++show_progress;
  }
  return true;
}

} /* namespace ComPWA */
