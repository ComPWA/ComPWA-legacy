#define BOOST_TEST_MODULE RootGeneratorTest

#include "Data/Root/RootGenerator.hpp"
#include <boost/test/unit_test.hpp>
#include <chrono>

BOOST_AUTO_TEST_SUITE(RootData)

unsigned int degreeOfDifferencePrecision(double difference,
                                         unsigned int maxdegree) {
  unsigned int index = difference / std::numeric_limits<double>::epsilon();
  if (index > maxdegree)
    index = maxdegree;
  return index;
}

void checkScenario(const ComPWA::FourMomentum &CMSP4,
                   const std::vector<double> &masses,
                   std::pair<double, double> EpsilonTolerances,
                   std::pair<double, double> TenEpsilonPercentages) {
  std::vector<int> FakePids;
  for (size_t i = 0; i < masses.size(); ++i) {
    FakePids.push_back(i);
  }
  auto EventGenerator =
      ComPWA::Data::Root::RootGenerator(CMSP4, masses, FakePids);

  double max_diff_masses(0.0);
  double max_diff_cms(0.0);
  std::vector<unsigned int> ToleranceDistributionCMS(30, 0);
  std::vector<unsigned int> ToleranceDistributionMasses(30, 0);
  unsigned int NumberOfEvents(100000);
  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(1234);

  auto EventCollect = EventGenerator.generate(NumberOfEvents, RandomGenerator);

  for (auto const &Event : EventCollect.Events) {
    auto const &FourVectors = Event.FourMomenta;

    for (unsigned int j = 0; j < masses.size(); ++j) {
      double tempdiff_mass(
          std::abs(FourVectors[j].invariantMass() - (double)masses[j]));
      if (tempdiff_mass > max_diff_masses) {
        // std::cout<<"new maxdiff!!!\n";
        max_diff_masses = tempdiff_mass;
      }
      unsigned int indexm = degreeOfDifferencePrecision(
          tempdiff_mass, ToleranceDistributionCMS.size() - 1);
      ++ToleranceDistributionMasses[indexm];
    }
    double tempdiff(
        std::abs(calculateInvariantMass(Event) - CMSP4.invariantMass()));
    if (tempdiff > max_diff_cms) {
      max_diff_cms = tempdiff;
    }
    unsigned int index = degreeOfDifferencePrecision(
        tempdiff, ToleranceDistributionCMS.size() - 1);
    ++ToleranceDistributionCMS[index];
  }

  BOOST_CHECK(max_diff_masses <
              EpsilonTolerances.first * std::numeric_limits<double>::epsilon());
  BOOST_CHECK(max_diff_cms < EpsilonTolerances.second *
                                 std::numeric_limits<double>::epsilon());

  // at 10*epsilon, we should have 95% of the events
  unsigned int count_masses(0);
  unsigned int count_cms(0);
  for (unsigned int i = 0; i < 10; ++i) {
    count_masses += ToleranceDistributionMasses[i];
    count_cms += ToleranceDistributionCMS[i];
  }

  BOOST_CHECK(100.0 * count_masses / (masses.size() * NumberOfEvents) >
              TenEpsilonPercentages.first);
  BOOST_CHECK(100.0 * count_cms / NumberOfEvents >
              TenEpsilonPercentages.second);
}

BOOST_AUTO_TEST_CASE(RootGeneratorPrecisionTest) {
  checkScenario(ComPWA::FourMomentum(0.0, 0.0, 0.0, 3.0), {0.2, 0.0},
                std::make_pair(1e9, 25), std::make_pair(85, 95));
  checkScenario(ComPWA::FourMomentum(0.0, 0.0, 0.0, 3.0), {0.2, 0.2, 0.0},
                std::make_pair(1e9, 25), std::make_pair(85, 95));
  checkScenario(ComPWA::FourMomentum(0.0, 0.0, 0.0, 3.0), {0.9, 0.9, 0.9},
                std::make_pair(25, 25), std::make_pair(95, 95));
  checkScenario(ComPWA::FourMomentum(0.0, 0.0, 0.0, 3.0), {2.0, 0.0, 0.0},
                std::make_pair(1e9, 25), std::make_pair(70, 95));
  checkScenario(ComPWA::FourMomentum(0.0, 0.0, 0.0, 3.0), {0.0, 0.0, 0.0},
                std::make_pair(1e9, 1e9), std::make_pair(60, 90));
  checkScenario(ComPWA::FourMomentum(0.0, 0.0, 0.0, 0.2), {0.0, 0.0, 0.0},
                std::make_pair(1e8, 1e8), std::make_pair(60, 95));
  checkScenario(ComPWA::FourMomentum(0.0, 0.0, 0.0, 4.0), {0.1, 0.5, 0.2, 0.3},
                std::make_pair(100, 25), std::make_pair(90, 95));
  checkScenario(ComPWA::FourMomentum(0.0, 0.0, 0.0, 4.0),
                {0.1, 0.5, 0.2, 0.3, 0.1}, std::make_pair(100, 25),
                std::make_pair(90, 95));
  checkScenario(ComPWA::FourMomentum(0.0, 0.0, 0.0, 4.0),
                {0.1, 0.5, 0.2, 0.3, 0.1, 0.2}, std::make_pair(100, 25),
                std::make_pair(90, 95));
};
BOOST_AUTO_TEST_SUITE_END()
