// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Define Boost test module
#define BOOST_TEST_MODULE Physics

#include "Physics/PhspVolume.hpp"
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace ComPWA::Physics;
using namespace std::chrono;

class Stopwatch {
public:
  Stopwatch() { start(); }
  void start() { t1 = steady_clock::now(); }
  double ns() {
    t2 = steady_clock::now();
    auto ns = (t2 - t1).count();
    t1 = steady_clock::now();
    return ns;
  }
  double ms() { return ns() / 1e6; }
  double s() { return ns() / 1e9; }

private:
  steady_clock::time_point t1;
  steady_clock::time_point t2;
};

double ComputeRelativeDifference(double value1, double value2) {
  double relDiff = value1 - value2;
  relDiff /= value1;
  relDiff *= 100.;
  relDiff = std::abs(relDiff);
  return relDiff;
}

class BenchmarkTable;

class Benchmark {
public:
  friend BenchmarkTable;
  Benchmark() {}
  Benchmark(std::string &TestDescription_, double ActualVolume_)
      : TestDescription(TestDescription_), ActualVolume(ActualVolume_) {}
  void SetActualVolume(double ActualVolume_) { ActualVolume = ActualVolume_; }
  void Compute(double ISMass_, std::vector<double> FSMasses_,
               size_t NumEvaluations_) {
    double s = ISMass_ * ISMass_;
    Stopwatch sw;
    auto vol = PhspVolume(s, FSMasses_, NumEvaluations_);
    CPUtime = sw.ms();
    Volume = vol.first;
    VolumeError = vol.second;
    RelativeDifference = ComputeRelativeDifference(Volume, ActualVolume);
    NSteps = NumEvaluations_;
  }
  void Print() {
    // Print description column
    std::cout << std::setw(wDescr) << std::left << TestDescription << " | ";
    // Print actual volume
    std::cout << std::setw(wActVol) << std::right << ActualVolume << " | ";
    // Print volume
    std::cout << std::setw(wVol) << std::right << Volume;
    // Print volume error
    if (VolumeError > 0.)
      std::cout << " +/- ";
    else
      std::cout << "     ";
    std::cout << std::setw(wVolErr) << std::left;
    if (VolumeError > 0.)
      std::cout << VolumeError;
    else
      std::cout << "";
    std::cout << " | ";
    // Print relative difference
    std::cout << std::setw(wRelDiff) << std::right;
    if (RelativeDifference > 0.)
      std::cout << RelativeDifference;
    else
      std::cout << "(equal)";
    std::cout << " | ";
    // Print number of iterations
    std::cout << std::setw(wNSteps) << std::right;
    if (NSteps > 1)
      std::cout << NSteps;
    else
      std::cout << "NA";
    std::cout << " | ";
    // Print CPU time
    std::cout << std::setw(wCPU) << std::right << CPUtime;
    std::cout << std::endl;
  }

private:
  std::string TestDescription;
  double ActualVolume;
  double Volume;
  double VolumeError;
  double RelativeDifference;
  double CPUtime;
  size_t NSteps;
  static const size_t wDescr = 12;
  static const size_t wActVol = 10;
  static const size_t wVol = 10;
  static const size_t wVolErr = 11;
  static const size_t wRelDiff = 12;
  static const size_t wNSteps = 10;
  static const size_t wCPU = 13;
};

class BenchmarkTable {
public:
  void AddBenchmark(std::string TestDescription_, double ISMass_,
                    std::vector<double> FSMasses_, size_t NumEvaluations_,
                    double ActualVolume_, double ToleranceInPercent) {
    Benchmark benchmark(TestDescription_, ActualVolume_);
    benchmark.Compute(ISMass_, FSMasses_, NumEvaluations_);
    BOOST_CHECK_CLOSE(benchmark.Volume, benchmark.ActualVolume,
                      ToleranceInPercent);
    Benchmarks.push_back(benchmark);
  }
  void Print() const {
    std::cout << std::endl;
    // Print header
    std::cout << std::left;
    std::cout << std::setw(Benchmark::wDescr) << "Test name"
              << " | ";
    std::cout << std::setw(Benchmark::wActVol) << "Actual vol"
              << " | ";
    std::cout << std::setw(Benchmark::wVol + Benchmark::wVolErr + 5) << "Volume"
              << " | ";
    std::cout << std::setw(Benchmark::wRelDiff) << "% difference"
              << " | ";
    std::cout << std::setw(Benchmark::wNSteps) << "N steps"
              << " | ";
    std::cout << std::setw(Benchmark::wCPU) << "CPU time (ms)";
    std::cout << std::endl;
    // Print line
    std::cout << std::setfill('-')
              << std::setw(Benchmark::wDescr + Benchmark::wActVol +
                           Benchmark::wVol + Benchmark::wVolErr +
                           Benchmark::wRelDiff + Benchmark::wNSteps +
                           Benchmark::wCPU + 20)
              << "" << std::endl;
    std::cout << std::setfill(' ');
    // Print rows
    for (auto bm : Benchmarks)
      bm.Print();
  }

private:
  std::vector<Benchmark> Benchmarks;
};

// Define Boost test suite
BOOST_AUTO_TEST_SUITE(Physics)

BOOST_AUTO_TEST_CASE(KallenFunction_test) {
  double RelativeTolerance(1e-6);
  // Test values
  BOOST_CHECK_CLOSE(KallenFunction(-1., 2., 3.), 12., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(-2., 3., 4.), 33., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(-3., 4., 5.), 64., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(-4., 5., 6.), 105., RelativeTolerance);
  // Test symmetry
  BOOST_CHECK_CLOSE(KallenFunction(2., 3., -1.), 12., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(3., -1., 2.), 12., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(3., 2., -1.), 12., RelativeTolerance);
}

BOOST_AUTO_TEST_CASE(PhspVolumeTest_DtoKKK) {
  // * Define initial state and final state masses
  double m0 = 1.86483;  // D0
  double m1 = 0.497611; // K0S
  double m2 = 0.493677; // K-
  double m3 = 0.493677; // K+
  std::vector<double> m1m2 = {m1, m2};
  std::vector<double> m1m2m3 = {m1, m2, m3};
  // * Test PhspVolume for 2 particles
  double actVol = 5.32194;
  BenchmarkTable table;
  table.AddBenchmark("D -> KK", m0, m1m2, 1, actVol, 1.);
  // * Test PhspVolume for 3 particles
  actVol = 3.0844;
  table.AddBenchmark("D -> KKK", m0, m1m2m3, 1e1, actVol, 5.);
  table.AddBenchmark("D -> KKK", m0, m1m2m3, 1e2, actVol, 1.);
  table.AddBenchmark("D -> KKK", m0, m1m2m3, 1e3, actVol, 1e-1);
  table.AddBenchmark("D -> KKK", m0, m1m2m3, 1e4, actVol, 1e-2);
  table.AddBenchmark("D -> KKK", m0, m1m2m3, 1e5, actVol, 1e-3);
  table.AddBenchmark("D -> KKK", m0, m1m2m3, 1e6, actVol, 1e-3);
  table.Print();
}

BOOST_AUTO_TEST_CASE(PhspVolumeTest_ngamma) {
  // * Define initial state and final state masses
  double m0 = 1.;
  double s = m0 * m0;
  std::vector<double> g2 = {0., 0.};
  std::vector<double> g3 = {0., 0., 0.};
  std::vector<double> g4 = {0., 0., 0., 0.};
  std::vector<double> g5 = {0., 0., 0., 0., 0.};
  // * Compute expected volumes
  double vol_2g = 2 * M_PI;
  double vol_3g = std::pow(M_PI, 2) * s;
  double vol_4g = std::pow(M_PI, 3) * std::pow(s, 2) / 6.;
  double vol_5g = std::pow(M_PI, 4) * std::pow(s, 3) / 72.;
  // * Test PhspVolume for e+e- --> n gamma
  BenchmarkTable table;
  table.AddBenchmark("ee -> 2g", m0, g2, 1, vol_2g, 1.);
  table.AddBenchmark("ee -> 3g", m0, g3, 1e4, vol_3g, 1.);
  table.AddBenchmark("ee -> 4g", m0, g4, 2e3, vol_4g, 1.);
  table.AddBenchmark("ee -> 5g", m0, g5, 300, vol_5g, 1.);
  table.Print();
}

BOOST_AUTO_TEST_SUITE_END()