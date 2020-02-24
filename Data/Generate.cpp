#include <algorithm>

#include "Core/Exceptions.hpp"
#include "Core/Generator.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Core/Random.hpp"
#include "Data/DataSet.hpp"
#include "Data/Generate.hpp"

#include "ThirdParty/parallelstl/include/pstl/algorithm"
#include "ThirdParty/parallelstl/include/pstl/execution"

namespace ComPWA {
namespace Data {

inline double uniform(double random, double min, double max) {
  return random * (max - min) + min;
}

std::tuple<std::vector<ComPWA::Event>, double>
generateBunch(unsigned int EventBunchSize, const ComPWA::Kinematics &Kinematics,
              ComPWA::Intensity &Intensity,
              ComPWA::UniformRealNumberGenerator &RandomGenerator,
              double generationMaxValue,
              std::vector<ComPWA::Event>::const_iterator PhspStartIterator,
              std::vector<ComPWA::Event>::const_iterator PhspTrueStartIterator,
              bool InverseIntensityWeighting = false) {

  std::vector<ComPWA::Event> SelectedEvents;

  auto TempDataSet = Kinematics.convert(std::vector<Event>(
      PhspTrueStartIterator, PhspTrueStartIterator + EventBunchSize));

  // evaluate the intensity
  auto Intensities = Intensity.evaluate(TempDataSet.Data);

  // multiply with event weights and set events outside of phase space (not
  // finite values = nan, inf) to zero
  std::vector<double> WeightedIntensities(Intensities.size());
  std::transform(pstl::execution::par_unseq, Intensities.begin(),
                 Intensities.end(), TempDataSet.Weights.begin(),
                 WeightedIntensities.begin(),
                 [](double intensity, double weight) {
                   if (std::isfinite(intensity)) {
                     return intensity * weight;
                   } else {
                     return 0.0;
                   }
                 });

  // determine maximum
  double BunchMax(*std::max_element(pstl::execution::par_unseq,
                                    WeightedIntensities.begin(),
                                    WeightedIntensities.end()));

  // restart generation if we got above the current maximum
  if (BunchMax > generationMaxValue) {
    return std::make_tuple(SelectedEvents, BunchMax);
  }

  // do hit and miss
  // first generate random numbers (no multithreading here, to ensure
  // deterministic behavior independent on the number of threads)
  std::vector<double> RandomNumbers;
  RandomNumbers.reserve(WeightedIntensities.size());
  std::generate_n(std::back_inserter(RandomNumbers), WeightedIntensities.size(),
                  [&RandomGenerator, generationMaxValue]() {
                    return uniform(RandomGenerator(), 0, generationMaxValue);
                  });

  if (InverseIntensityWeighting) {
    for (unsigned int i = 0; i < WeightedIntensities.size(); ++i) {
      if (RandomNumbers[i] < WeightedIntensities[i]) {
        SelectedEvents.push_back(ComPWA::Event{PhspStartIterator->FourMomenta,
                                               1.0 / Intensities[i]});
      }
      ++PhspStartIterator;
    }
  } else {
    for (unsigned int i = 0; i < WeightedIntensities.size(); ++i) {
      if (RandomNumbers[i] < WeightedIntensities[i]) {
        SelectedEvents.push_back(
            ComPWA::Event{PhspStartIterator->FourMomenta, 1.0});
      }
      ++PhspStartIterator;
    }
  }

  return std::make_tuple(SelectedEvents, BunchMax);
}

ComPWA::EventList
generate(unsigned int NumberOfEvents, const ComPWA::Kinematics &Kinematics,
         const ComPWA::PhaseSpaceEventGenerator &Generator,
         ComPWA::Intensity &Intensity,
         ComPWA::UniformRealNumberGenerator &RandomGenerator) {
  if (NumberOfEvents <= 0)
    return ComPWA::EventList{};

  std::vector<Event> Events;
  // initialize generator output vector
  unsigned int EventBunchSize(5000);
  Events.reserve(NumberOfEvents);
  unsigned int TotalGeneratedEvents(0);

  double SafetyMargin(0.05);
  double generationMaxValue(0.0);
  unsigned int initialSeed = RandomGenerator.getSeed();

  LOG(INFO) << "Generating hit-and-miss sample: [" << NumberOfEvents
            << " events] ";
  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    std::vector<ComPWA::Event> TempEvents;
    TotalGeneratedEvents += EventBunchSize;
    // generate events
    std::generate_n(std::back_inserter(TempEvents), EventBunchSize,
                    [&Generator, &RandomGenerator]() {
                      return Generator.generate(RandomGenerator);
                    });

    auto Bunch = generateBunch(EventBunchSize, Kinematics, Intensity,
                               RandomGenerator, generationMaxValue,
                               TempEvents.begin(), TempEvents.begin());

    std::vector<Event> BunchEvents = std::get<0>(Bunch);
    double MaximumWeight = std::get<1>(Bunch);

    // restart generation if we got above the current maximum
    if (MaximumWeight > generationMaxValue) {
      generationMaxValue = (1.0 + SafetyMargin) * MaximumWeight;

      if (Events.size() > 0) {
        Events.clear();
        RandomGenerator.setSeed(initialSeed);
        bar = ComPWA::ProgressBar(NumberOfEvents);
        LOG(INFO) << "Tools::generate() | Error in HitMiss "
                     "procedure: Maximum value of random number generation "
                     "smaller then amplitude maximum! We raise the maximum "
                     "to "
                  << generationMaxValue << " value and restart generation!";
        continue;
      }
    }

    size_t AmountToAppend(BunchEvents.size());
    if (Events.size() + BunchEvents.size() > NumberOfEvents) {
      AmountToAppend = NumberOfEvents - Events.size();
    }

    Events.insert(
        Events.end(), std::make_move_iterator(BunchEvents.begin()),
        std::make_move_iterator(BunchEvents.begin() + AmountToAppend));
    bar.next(AmountToAppend);

    if (Events.size() == NumberOfEvents)
      break;
  }
  LOG(INFO) << "Successfully generated " << NumberOfEvents
            << " with an efficiency of "
            << 1.0 * NumberOfEvents / TotalGeneratedEvents;

  return EventList{Kinematics.getFinalStatePIDs(), Events};
}

EventList generate(unsigned int NumberOfEvents,
                   const ComPWA::Kinematics &Kinematics,
                   ComPWA::UniformRealNumberGenerator &RandomGenerator,
                   ComPWA::Intensity &Intensity,
                   const std::vector<ComPWA::Event> &phsp,
                   const std::vector<ComPWA::Event> &phspTrue) {
  // Doing some checks
  if (NumberOfEvents <= 0)
    throw std::runtime_error("Tools::generate() negative number of events: " +
                             std::to_string(NumberOfEvents));

  if (phspTrue.size() != phsp.size())
    throw std::runtime_error(
        "Tools::generate() | We have a sample of true "
        "phsp events, but the sample size doesn't match that one of "
        "the phsp sample!");

  std::vector<ComPWA::Event> events;
  events.reserve(NumberOfEvents);

  double SafetyMargin(0.05);

  double maxSampleWeight(ComPWA::getMaximumSampleWeight(phsp));
  if (phspTrue.size()) {
    double temp_maxweight(ComPWA::getMaximumSampleWeight(phspTrue));
    if (temp_maxweight > maxSampleWeight)
      maxSampleWeight = temp_maxweight;
  }

  if (maxSampleWeight <= 0.0)
    throw std::runtime_error("Tools::generate() Sample maximum value is zero!");
  double generationMaxValue(maxSampleWeight * (1.0 + SafetyMargin));
  unsigned int initialSeed = RandomGenerator.getSeed();

  LOG(INFO) << "Tools::generate() | Using " << generationMaxValue
            << " as maximum value of the intensity.";

  auto const &PhspEvents = phsp;
  unsigned int LastIndex(PhspEvents.size() - 1);

  unsigned int EventBunchSize(5000);
  if (PhspEvents.size() < EventBunchSize)
    EventBunchSize = PhspEvents.size();

  auto CurrentStartIterator = PhspEvents.begin();
  auto CurrentTrueStartIterator = PhspEvents.begin();
  if (phspTrue.size())
    CurrentTrueStartIterator = phspTrue.begin();
  unsigned int CurrentStartIndex(0);
  LOG(INFO) << "Generating hit-and-miss sample: [" << NumberOfEvents
            << " events] ";
  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    if (CurrentStartIndex + EventBunchSize > LastIndex)
      EventBunchSize = LastIndex - CurrentStartIndex;

    auto Bunch = generateBunch(EventBunchSize, Kinematics, Intensity,
                               RandomGenerator, generationMaxValue,
                               CurrentStartIterator, CurrentTrueStartIterator);

    std::vector<Event> BunchEvents = std::get<0>(Bunch);
    double MaximumWeight = std::get<1>(Bunch);

    // restart generation if we got above the current maximum
    if (MaximumWeight > generationMaxValue) {
      generationMaxValue = (1.0 + SafetyMargin) * MaximumWeight;
      LOG(INFO) << "We raise the maximum to " << generationMaxValue;
      if (events.size() > 0) {
        events.clear();
        RandomGenerator.setSeed(initialSeed);
        CurrentStartIterator = PhspEvents.begin();
        CurrentTrueStartIterator = PhspEvents.begin();
        if (phspTrue.size())
          CurrentTrueStartIterator = phspTrue.begin();
        CurrentStartIndex = 0;
        bar = ComPWA::ProgressBar(NumberOfEvents);
        LOG(INFO) << "Tools::generate() | Error in HitMiss "
                     "procedure: Maximum value of random number generation "
                     "smaller then amplitude maximum! Restarting generation!";
      }
      continue;
    }

    size_t AmountToAppend(BunchEvents.size());
    if (events.size() + BunchEvents.size() > NumberOfEvents) {
      AmountToAppend = NumberOfEvents - events.size();
    }
    events.insert(
        events.end(), std::make_move_iterator(BunchEvents.begin()),
        std::make_move_iterator(BunchEvents.begin() + AmountToAppend));

    bar.next(AmountToAppend);

    if (events.size() == NumberOfEvents)
      break;

    // increment true iterator
    std::advance(CurrentTrueStartIterator, EventBunchSize);
    std::advance(CurrentStartIterator, EventBunchSize);
    CurrentStartIndex += EventBunchSize;

    if (CurrentStartIndex >= LastIndex)
      break;
  }
  double gen_eff = (double)events.size() / NumberOfEvents;
  if (CurrentStartIndex > NumberOfEvents) {
    gen_eff = (double)events.size() / CurrentStartIndex;
  }
  LOG(INFO) << "Efficiency of toy MC generation: " << gen_eff;

  return EventList{Kinematics.getFinalStatePIDs(), events};
}

std::vector<Event>
generatePhsp(unsigned int nEvents,
             const ComPWA::PhaseSpaceEventGenerator &Generator,
             ComPWA::UniformRealNumberGenerator &RandomGenerator) {
  std::vector<Event> Events;

  LOG(INFO) << "Generating phase-space MC: [" << nEvents << " events] ";

  ComPWA::ProgressBar bar(nEvents);
  for (unsigned int i = 0; i < nEvents; ++i) {
    ComPWA::Event tmp = Generator.generate(RandomGenerator);
    double ampRnd = RandomGenerator();
    if (ampRnd > tmp.Weight) {
      --i;
      continue;
    }

    // Reset weights: weights are taken into account by hit&miss. The
    // resulting sample is therefore unweighted
    Events.push_back(ComPWA::Event{tmp.FourMomenta, 1.0});
    bar.next();
  }
  return Events;
}

EventList generateImportanceSampledPhsp(
    unsigned int NumberOfEvents, const ComPWA::Kinematics &Kinematics,
    const ComPWA::PhaseSpaceEventGenerator &Generator,
    ComPWA::Intensity &Intensity,
    ComPWA::UniformRealNumberGenerator &RandomGenerator) {
  std::vector<ComPWA::Event> Events;
  if (NumberOfEvents <= 0)
    return EventList{};
  // initialize generator output vector
  unsigned int EventBunchSize(5000);
  Events.reserve(NumberOfEvents);

  double SafetyMargin(0.05);
  double generationMaxValue(0.0);
  unsigned int initialSeed = RandomGenerator.getSeed();

  LOG(INFO)
      << "Generating phase space sample (hit-and-miss importance sampled): ["
      << NumberOfEvents << " events] ";
  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    std::vector<ComPWA::Event> TempEvents;
    // generate events
    std::generate_n(std::back_inserter(TempEvents), EventBunchSize,
                    [&Generator, &RandomGenerator]() {
                      return Generator.generate(RandomGenerator);
                    });

    std::vector<Event> BunchEvents;
    double MaximumWeight;

    std::tie(BunchEvents, MaximumWeight) = generateBunch(
        EventBunchSize, Kinematics, Intensity, RandomGenerator,
        generationMaxValue, TempEvents.begin(), TempEvents.begin(), true);

    // restart generation if we got above the current maximum
    if (MaximumWeight > generationMaxValue) {
      generationMaxValue = (1.0 + SafetyMargin) * MaximumWeight;
      if (Events.size() > 0) {
        Events.clear();
        RandomGenerator.setSeed(initialSeed);
        bar = ComPWA::ProgressBar(NumberOfEvents);
        LOG(INFO)
            << "Tools::generateImportanceSampledPhsp() | Error in HitMiss "
               "procedure: Maximum value of random number generation "
               "smaller then amplitude maximum! We raise the maximum "
               "to "
            << generationMaxValue << " value and restart generation!";
        continue;
      }
    }

    size_t AmountToAppend(BunchEvents.size());
    if (Events.size() + BunchEvents.size() > NumberOfEvents) {
      AmountToAppend = NumberOfEvents - Events.size();
    }

    Events.insert(Events.end(), BunchEvents.begin(),
                  BunchEvents.begin() + AmountToAppend);

    bar.next(AmountToAppend);

    if (Events.size() == NumberOfEvents)
      break;
  }
  // replace with std::reduce once standard is moved to c++17
  double WeightSum(0.0);
  for (auto const &x : Events) {
    WeightSum += x.Weight;
  }

  // now just rescale the event weights so that sum(event weights) = # events
  double rescale_factor(NumberOfEvents / WeightSum);
  for (auto &evt : Events) {
    evt.Weight *= rescale_factor;
  }

  return EventList{Kinematics.getFinalStatePIDs(), Events};
}

} // namespace Data
} // namespace ComPWA
