// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Some useful function for Monte-Carlo generation.
///

#ifndef Generate_hpp
#define Generate_hpp

#include "Core/ProgressBar.hpp"
#include "Core/Generator.hpp"
#include "DataReader/Data.hpp"
#include "Core/AmpIntensity.hpp"
#include "Core/Kinematics.hpp"
#include "Tools/Integration.hpp"

namespace ComPWA {
namespace Tools {

inline bool Generate(int number, std::shared_ptr<ComPWA::Kinematics> kin,
                     std::shared_ptr<ComPWA::Generator> gen,
                     std::shared_ptr<ComPWA::AmpIntensity> amp,
                     std::shared_ptr<ComPWA::DataReader::Data> data,
                     std::shared_ptr<ComPWA::DataReader::Data> phsp,
                     std::shared_ptr<ComPWA::DataReader::Data> phspTrue =
                         std::shared_ptr<ComPWA::DataReader::Data>()) {

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
  double generationMaxValue = 5 * maxSampleWeight;

  unsigned int initialSeed = gen->GetSeed();
  unsigned int totalCalls = 0;

  unsigned int limit = 100000000; // set large limit, should never be reached;
  if (phsp) {
    limit = phsp->GetNEvents();
    generationMaxValue *= ComPWA::Tools::Maximum(kin, amp, phsp);
  }

  LOG(trace) << "RunMananger::gen() | Using " << generationMaxValue
             << " as maximum value of the intensity.";

  ComPWA::Event evt;     // event that we fill into generated sample
  ComPWA::Event evtTrue; // event that is used to evalutate amplitude
  ComPWA::progressBar bar(number);
  if (number <= 0)
    bar = ComPWA::progressBar(limit);
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
    ComPWA::dataPoint point;
    try {
      kin->EventToDataPoint(evtTrue, point);
    } catch (ComPWA::BeyondPhsp &ex) { // event outside phase, remove
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
      bar = ComPWA::progressBar(number);
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
  if (data->GetNEvents() < number) {
    std::cout << std::endl;
    LOG(error) << "RunManager::gen() | Not able to generate " << number
               << " events. Phsp sample too small. Current size "
                  "of sample is now "
               << data->GetNEvents();
  }

  if (!totalCalls)
    throw std::runtime_error("RunManager::gen() | Number of calls is zero! "
                             "There ust be something wrong!");

  double genEff = (double)data->GetNEvents() / totalCalls;
  LOG(info) << "Efficiency of toy MC generation: " << genEff << ".";

  return true;
}

inline bool GeneratePhsp(int nEvents, std::shared_ptr<ComPWA::Generator> gen,
                         std::shared_ptr<ComPWA::DataReader::Data> sample) {
  if (nEvents == 0)
    return 0;
  if (!sample)
    throw std::runtime_error("RunManager: generatePhsp() | "
                             "No phase-space sample set");
  if (sample->GetNEvents() > 0)
    throw std::runtime_error("RunManager: generatePhsp() | "
                             "Dataset not empty! abort!");

  LOG(info) << "Generating phase-space MC: [" << nEvents << " events] ";

  ComPWA::progressBar bar(nEvents);
  for (unsigned int i = 0; i < nEvents; i++) {
    if (i > 0)
      i--;
    ComPWA::Event tmp;
    gen->Generate(tmp);
    double ampRnd = gen->GetUniform(0, 1);
    if (ampRnd > tmp.GetWeight())
      continue;

    // Reset weights: weights are taken into account by hit&miss. The
    // resulting sample is therefore unweighted
    tmp.SetWeight(1.);

    tmp.SetEfficiency(1.);
    i++;
    sample->PushEvent(tmp); // unfortunatly not thread safe
    bar.nextEvent();
  }
  return true;
}

} // ns::Tools
} // ns::ComPWA
#endif
