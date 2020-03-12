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

std::tuple<EventCollection, double>
generateBunch(unsigned int EventBunchSize, const ComPWA::Kinematics &Kinematics,
              ComPWA::Intensity &Intensity,
              ComPWA::UniformRealNumberGenerator &RandomGenerator,
              double generationMaxValue,
              std::vector<ComPWA::Event>::const_iterator PhspStartIterator,
              std::vector<ComPWA::Event>::const_iterator PhspTrueStartIterator,
              bool InverseIntensityWeighting = false) {

  EventCollection SelectedEvents{Kinematics.getFinalStatePIDs()};

  auto NewEvents = EventCollection{
      Kinematics.getFinalStatePIDs(),
      {PhspTrueStartIterator, PhspTrueStartIterator + EventBunchSize}};
  auto TempDataSet = Kinematics.convert(NewEvents);

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
        SelectedEvents.Events.push_back(
            {PhspStartIterator->FourMomenta, 1.0 / Intensities[i]});
      }
      ++PhspStartIterator;
    }
  } else {
    for (unsigned int i = 0; i < WeightedIntensities.size(); ++i) {
      if (RandomNumbers[i] < WeightedIntensities[i]) {
        SelectedEvents.Events.push_back({PhspStartIterator->FourMomenta, 1.0});
      }
      ++PhspStartIterator;
    }
  }

  return std::make_tuple(SelectedEvents, BunchMax);
}

EventCollection generate(unsigned int NumberOfEvents,
                         const ComPWA::Kinematics &Kinematics,
                         const ComPWA::PhaseSpaceEventGenerator &Generator,
                         ComPWA::Intensity &Intensity,
                         ComPWA::UniformRealNumberGenerator &RandomGenerator) {
  EventCollection GeneratedEvents{Kinematics.getFinalStatePIDs()};
  if (NumberOfEvents <= 0)
    return GeneratedEvents;

  // initialize generator output vector
  unsigned int EventBunchSize(5000);
  GeneratedEvents.Events.reserve(NumberOfEvents);
  unsigned int TotalGeneratedEvents(0);

  double SafetyMargin(0.05);
  double generationMaxValue(0.0);
  unsigned int initialSeed = RandomGenerator.getSeed();

  LOG(INFO) << "Generating hit-and-miss sample: [" << NumberOfEvents
            << " events] ";
  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    EventCollection TempEvents =
        Generator.generate(EventBunchSize, RandomGenerator);

    TotalGeneratedEvents += EventBunchSize;

    auto Bunch =
        generateBunch(EventBunchSize, Kinematics, Intensity, RandomGenerator,
                      generationMaxValue, TempEvents.Events.begin(),
                      TempEvents.Events.begin());

    EventCollection BunchEvents = std::get<0>(Bunch);
    double MaximumWeight = std::get<1>(Bunch);

    // restart generation if we got above the current maximum
    if (MaximumWeight > generationMaxValue) {
      generationMaxValue = (1.0 + SafetyMargin) * MaximumWeight;

      if (GeneratedEvents.Events.size() > 0) {
        GeneratedEvents.Events.clear();
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

    size_t AmountToAppend(BunchEvents.Events.size());
    if (GeneratedEvents.Events.size() + BunchEvents.Events.size() >
        NumberOfEvents) {
      AmountToAppend = NumberOfEvents - GeneratedEvents.Events.size();
    }

    GeneratedEvents.Events.insert(
        GeneratedEvents.Events.end(),
        std::make_move_iterator(BunchEvents.Events.begin()),
        std::make_move_iterator(BunchEvents.Events.begin() + AmountToAppend));
    bar.next(AmountToAppend);

    if (GeneratedEvents.Events.size() == NumberOfEvents)
      break;
  }
  LOG(INFO) << "Successfully generated " << NumberOfEvents
            << " with an efficiency of "
            << 1.0 * NumberOfEvents / TotalGeneratedEvents;

  return GeneratedEvents;
}

EventCollection generate(unsigned int NumberOfEvents,
                         const ComPWA::Kinematics &Kinematics,
                         ComPWA::UniformRealNumberGenerator &RandomGenerator,
                         ComPWA::Intensity &Intensity,
                         const EventCollection &PhspSample,
                         const EventCollection &PhspSampleTrue) {
  // Doing some checks
  if (NumberOfEvents <= 0)
    throw std::runtime_error("Tools::generate() negative number of events: " +
                             std::to_string(NumberOfEvents));

  if (PhspSampleTrue.Events.size() != PhspSample.Events.size())
    throw std::runtime_error(
        "Tools::generate() | We have a sample of true "
        "phsp events, but the sample size doesn't match that one of "
        "the phsp sample!");

  EventCollection GeneratedEvents{Kinematics.getFinalStatePIDs()};
  GeneratedEvents.Events.resize(NumberOfEvents);

  double SafetyMargin(0.05);

  double maxSampleWeight(ComPWA::getMaximumSampleWeight(PhspSample));
  if (PhspSampleTrue.Events.size()) {
    double temp_maxweight(ComPWA::getMaximumSampleWeight(PhspSampleTrue));
    if (temp_maxweight > maxSampleWeight)
      maxSampleWeight = temp_maxweight;
  }

  if (maxSampleWeight <= 0.0)
    throw std::runtime_error("Tools::generate() Sample maximum value is zero!");
  double generationMaxValue(maxSampleWeight * (1.0 + SafetyMargin));
  unsigned int initialSeed = RandomGenerator.getSeed();

  LOG(INFO) << "Tools::generate() | Using " << generationMaxValue
            << " as maximum value of the intensity.";

  auto const &PhspEvents = PhspSample;
  unsigned int LastIndex(PhspEvents.Events.size() - 1);

  unsigned int EventBunchSize(5000);
  if (PhspEvents.Events.size() < EventBunchSize)
    EventBunchSize = PhspEvents.Events.size();

  auto CurrentStartIterator = PhspEvents.Events.begin();
  auto CurrentTrueStartIterator = PhspEvents.Events.begin();
  if (PhspSampleTrue.Events.size())
    CurrentTrueStartIterator = PhspSampleTrue.Events.begin();
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

    EventCollection BunchEvents = std::get<0>(Bunch);
    double MaximumWeight = std::get<1>(Bunch);

    // restart generation if we got above the current maximum
    if (MaximumWeight > generationMaxValue) {
      generationMaxValue = (1.0 + SafetyMargin) * MaximumWeight;
      LOG(INFO) << "We raise the maximum to " << generationMaxValue;
      if (GeneratedEvents.Events.size() > 0) {
        GeneratedEvents.Events.clear();
        RandomGenerator.setSeed(initialSeed);
        CurrentStartIterator = PhspEvents.Events.begin();
        CurrentTrueStartIterator = PhspEvents.Events.begin();
        if (PhspSampleTrue.Events.size())
          CurrentTrueStartIterator = PhspSampleTrue.Events.begin();
        CurrentStartIndex = 0;
        bar = ComPWA::ProgressBar(NumberOfEvents);
        LOG(INFO) << "Tools::generate() | Error in HitMiss "
                     "procedure: Maximum value of random number generation "
                     "smaller then amplitude maximum! Restarting generation!";
      }
      continue;
    }

    size_t AmountToAppend(BunchEvents.Events.size());
    if (GeneratedEvents.Events.size() + BunchEvents.Events.size() >
        NumberOfEvents) {
      AmountToAppend = NumberOfEvents - GeneratedEvents.Events.size();
    }
    GeneratedEvents.Events.insert(
        GeneratedEvents.Events.end(),
        std::make_move_iterator(BunchEvents.Events.begin()),
        std::make_move_iterator(BunchEvents.Events.begin() + AmountToAppend));

    bar.next(AmountToAppend);

    if (GeneratedEvents.Events.size() == NumberOfEvents)
      break;

    // increment true iterator
    std::advance(CurrentTrueStartIterator, EventBunchSize);
    std::advance(CurrentStartIterator, EventBunchSize);
    CurrentStartIndex += EventBunchSize;

    if (CurrentStartIndex >= LastIndex)
      break;
  }
  double gen_eff = (double)GeneratedEvents.Events.size() / NumberOfEvents;
  if (CurrentStartIndex > NumberOfEvents) {
    gen_eff = (double)GeneratedEvents.Events.size() / CurrentStartIndex;
  }
  LOG(INFO) << "Efficiency of toy MC generation: " << gen_eff;

  return GeneratedEvents;
}

EventCollection
generatePhsp(unsigned int NumberOfEvents,
             const ComPWA::PhaseSpaceEventGenerator &Generator,
             ComPWA::UniformRealNumberGenerator &RandomGenerator) {
  EventCollection GeneratedPhsp = Generator.generate(0, RandomGenerator);
  unsigned int EventBunchSize(5000);

  LOG(INFO) << "Generating phase-space MC: [" << NumberOfEvents << " events] ";

  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    EventCollection TmpEvents =
        Generator.generate(EventBunchSize, RandomGenerator);
    std::vector<double> RandomNumbers;
    RandomNumbers.reserve(EventBunchSize);
    std::generate_n(std::back_inserter(RandomNumbers), EventBunchSize,
                    [&RandomGenerator]() { return RandomGenerator(); });
    for (size_t i = 0; i < RandomNumbers.size(); ++i) {
      if (GeneratedPhsp.Events.size() == NumberOfEvents)
        break;
      if (RandomNumbers[i] > TmpEvents.Events[i].Weight) {
        continue;
      }

      // Reset weights: weights are taken into account by hit&miss. The
      // resulting sample is therefore unweighted
      GeneratedPhsp.Events.push_back(
          ComPWA::Event{TmpEvents.Events[i].FourMomenta, 1.0});
      bar.next();
    }
    if (GeneratedPhsp.Events.size() == NumberOfEvents)
      break;
  }
  return GeneratedPhsp;
}

EventCollection generateImportanceSampledPhsp(
    unsigned int NumberOfEvents, const ComPWA::Kinematics &Kinematics,
    const ComPWA::PhaseSpaceEventGenerator &Generator,
    ComPWA::Intensity &Intensity,
    ComPWA::UniformRealNumberGenerator &RandomGenerator) {
  EventCollection GeneratedEventList{Kinematics.getFinalStatePIDs()};
  if (NumberOfEvents <= 0)
    return GeneratedEventList;
  // initialize generator output vector
  unsigned int EventBunchSize(5000);
  GeneratedEventList.Events.reserve(NumberOfEvents);

  double SafetyMargin(0.05);
  double generationMaxValue(0.0);
  unsigned int initialSeed = RandomGenerator.getSeed();

  LOG(INFO)
      << "Generating phase space sample (hit-and-miss importance sampled): ["
      << NumberOfEvents << " events] ";
  ComPWA::ProgressBar bar(NumberOfEvents);
  while (true) {
    EventCollection TempEventList =
        Generator.generate(EventBunchSize, RandomGenerator);

    EventCollection BunchEvents{Kinematics.getFinalStatePIDs()};
    double MaximumWeight;

    std::tie(BunchEvents, MaximumWeight) =
        generateBunch(EventBunchSize, Kinematics, Intensity, RandomGenerator,
                      generationMaxValue, TempEventList.Events.begin(),
                      TempEventList.Events.begin(), true);

    // restart generation if we got above the current maximum
    if (MaximumWeight > generationMaxValue) {
      generationMaxValue = (1.0 + SafetyMargin) * MaximumWeight;
      if (GeneratedEventList.Events.size() > 0) {
        GeneratedEventList.Events.clear();
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

    size_t AmountToAppend(BunchEvents.Events.size());
    if (GeneratedEventList.Events.size() + BunchEvents.Events.size() >
        NumberOfEvents) {
      AmountToAppend = NumberOfEvents - GeneratedEventList.Events.size();
    }

    GeneratedEventList.Events.insert(
        GeneratedEventList.Events.end(), BunchEvents.Events.begin(),
        BunchEvents.Events.begin() + AmountToAppend);

    bar.next(AmountToAppend);

    if (GeneratedEventList.Events.size() == NumberOfEvents)
      break;
  }
  // replace with std::reduce once standard is moved to c++17
  double WeightSum(0.0);
  for (auto const &Event : GeneratedEventList.Events) {
    WeightSum += Event.Weight;
  }

  // now just rescale the event weights so that sum(event weights) = # events
  double rescale_factor(NumberOfEvents / WeightSum);
  for (auto &Event : GeneratedEventList.Events) {
    Event.Weight *= rescale_factor;
  }

  return GeneratedEventList;
}

} // namespace Data
} // namespace ComPWA
