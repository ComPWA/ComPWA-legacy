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

inline bool generate(int number, std::shared_ptr<ComPWA::Kinematics> kin,
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
    throw std::runtime_error("Tools::generate() negative number of events: " +
                             std::to_string((long double)number) +
                             ". And no phsp sample given!");
  if (!amp)
    throw std::runtime_error("Tools::generate() | Amplitude not valid");
  if (!gen)
    throw std::runtime_error("Tools::generate() | Generator not valid");
  if (!data)
    throw std::runtime_error("Tools::generate() | Sample not valid");
  if (data->numEvents() > 0)
    throw std::runtime_error("Tools::generate() | Sample not empty!");
  if (phspTrue && !phsp)
    throw std::runtime_error("Tools::generate() | We have a sample of true"
                             " phsp events, but no phsp sample!");
  if (phspTrue && phspTrue->numEvents() != phsp->numEvents())
    throw std::runtime_error(
        "Tools::generate() | We have a sample of true "
        "phsp events, but the sample size doesn't match that one of "
        "the phsp sample!");

  double maxSampleWeight = 1.0;
  if (phsp)
    maxSampleWeight = phsp->maximumWeight();
  if (phspTrue && phspTrue->maximumWeight() > maxSampleWeight)
    maxSampleWeight = phspTrue->maximumWeight();

  // Maximum value for random number generation. We introduce an arbitrary
  // factor of 5 to make sure that the maximum value is never reached.
  double generationMaxValue = 5 * maxSampleWeight;

  unsigned int initialSeed = gen->seed();
  unsigned int totalCalls = 0;

  unsigned int limit = 100000000; // set large limit, should never be reached;
  if (phsp) {
    limit = phsp->numEvents();
    generationMaxValue *= ComPWA::Tools::Maximum(kin, amp, phsp);
  }

  LOG(TRACE) << "Tools::generate() | Using " << generationMaxValue
             << " as maximum value of the intensity.";

  ComPWA::Event evt;     // event that we fill into generated sample
  ComPWA::Event evtTrue; // event that is used to evalutate amplitude
  ComPWA::ProgressBar bar(number);
  if (number <= 0)
    bar = ComPWA::ProgressBar(limit);
  for (unsigned int i = 0; i < limit; i++) {
    if (phsp && phspTrue) { // phsp and true sample is set
      evtTrue = phspTrue->event(i);
      evt = phsp->event(i);
    } else if (phsp) { // phsp sample is set
      evt = phsp->event(i);
      evtTrue = evt;
    } else { // otherwise generate event
      gen->generate(evt);
      evtTrue = evt;
    }
    if (number <= 0)
      bar.next();

    // use reconstructed position for weights
    double weight = evt.weight();

    // use true position for amplitude value
    ComPWA::DataPoint point;
    try {
      kin->convert(evtTrue, point);
    } catch (ComPWA::BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }

    totalCalls++;
    double ampRnd = gen->uniform(0, generationMaxValue);
    double AMPpdf = amp->intensity(point);

    // If maximum of intensity is reached we have to restart the procedure
    if (generationMaxValue < (AMPpdf * weight)) {
      i = 0;
      bar = ComPWA::ProgressBar(number);
      gen->setSeed(initialSeed);
      generationMaxValue = 2 * (AMPpdf * weight);
      data->clear();
      totalCalls = 0;
      LOG(TRACE) << "Tools::generate() | Error in HitMiss "
                    "procedure: Maximum value of random number generation "
                    "smaller then amplitude maximum! We raise the maximum "
                    "to "
                 << generationMaxValue << " value and restart generation!";
      continue;
    }

    if (ampRnd > (weight * AMPpdf))
      continue;

    // Fill event to sample
    // reset weights: the weights are taken into account by hit and miss. The
    // resulting sample is therefore unweighted
    evt.setWeight(1.);     // reset weight
    evt.setEfficiency(1.); // reset weight
    data->add(evt);

    if (number > 0)
      bar.next();
    
    // break if we have a sufficienct number of events
    if (data->numEvents() >= number)
      i = limit;
  }
  
  if (data->numEvents() < number) {
    std::cout << std::endl;
    LOG(ERROR) << "Tools::generate() | Not able to generate " << number
               << " events. Phsp sample too small. Current size "
                  "of sample is now "
               << data->numEvents();
  }

  if (!totalCalls)
    throw std::runtime_error("Tools::generate() | Number of calls is zero! "
                             "There must be something wrong!");

  double genEff = (double)data->numEvents() / totalCalls;
  LOG(INFO) << "Efficiency of toy MC generation: " << genEff << ".";

  return true;
}

inline bool generatePhsp(int nEvents, std::shared_ptr<ComPWA::Generator> gen,
                         std::shared_ptr<ComPWA::DataReader::Data> sample) {
  if (nEvents == 0)
    return 0;
  if (!sample)
    throw std::runtime_error("Tools::GeneratePhsp() | "
                             "No phase-space sample set");
  if (sample->numEvents() > 0)
    throw std::runtime_error("Tools::GeneratePhsp() | "
                             "Dataset not empty! abort!");

  LOG(INFO) << "Generating phase-space MC: [" << nEvents << " events] ";

  ComPWA::ProgressBar bar(nEvents);
  for (unsigned int i = 0; i < nEvents; i++) {
    if (i > 0)
      i--;
    ComPWA::Event tmp;
    gen->generate(tmp);
    double ampRnd = gen->uniform(0, 1);
    if (ampRnd > tmp.weight())
      continue;

    // Reset weights: weights are taken into account by hit&miss. The
    // resulting sample is therefore unweighted
    tmp.setWeight(1.);

    tmp.setEfficiency(1.);
    i++;
    sample->add(tmp); // unfortunatly not thread safe
    bar.next();
  }
  return true;
}

} // ns::Tools
} // ns::ComPWA
#endif
