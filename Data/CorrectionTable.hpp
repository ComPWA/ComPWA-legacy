// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Correction table.
///

#ifndef COMPWA_DATA_MOMENTUMSYSTEMATICS_HPP_
#define COMPWA_DATA_MOMENTUMSYSTEMATICS_HPP_

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ComPWA {
namespace Data {
///
/// \class CorrectionTable
/// Corrections in FourMomenta taken from a look-up-table and a weight is
/// calculated.
///
class CorrectionTable {
public:
  CorrectionTable(std::string t = "") : title(t) {}

  ~CorrectionTable() {}

  /// Get Data/MC difference
  /// Choose a charge for particle(=1) or anti-particle(=-1). A charge of 0
  /// averages over both values.
  double GetValue(int charge, double momentum) const;

  /// Get error of Data/MC difference
  /// Choose a charge for particle(=1) or anti-particle(=-1). A charge of 0
  /// averages over both values.
  double GetError(int charge, double momentum) const;

  /// Get MC/Data difference
  /// Choose a charge for particle(=1) or anti-particle(=-1). A charge of 0
  /// averages over both values.
  double GetInvValue(int charge, double momentum) const;

  /// Get error of MC/Data difference
  /// Choose a charge for particle(=1) or anti-particle(=-1). A charge of 0
  /// averages over both values.
  double GetInvError(int charge, double momentum) const;

  /// Add momentum bin with efficiency correction
  /// The efficiency is given as (epsilon_data/epsilon_mc - 1)
  /// \param binMin lower bin boundary
  /// \param binMax higher bin boundary
  /// \param s efficiency difference for particles
  /// \param sError uncertainty of efficiency difference for particles
  /// \param antiS efficiency difference for anti-particles
  /// \param antiSerror uncertainty of efficiency difference for anti-particles
  void Add(double binMin, double binMax, double s, double sError,
           double antiS = -999, double antiSerror = -999);

  /// Add momentum bin with efficiency correction
  /// The efficiency is given as (epsilon_mc/epsilon_data - 1)
  /// \param binMin lower bin boundary
  /// \param binMax higher bin boundary
  /// \param s efficiency difference for particles
  /// \param sError uncertainty of efficiency difference for particles
  /// \param antiS efficiency difference for anti-particles
  /// \param antiSerror uncertainty of efficiency difference for anti-particles
  void AddInv(double binMin, double binMax, double s, double sError,
              double antiS = -999, double antiSerror = -999);

  std::vector<std::pair<double, double>> bins() { return Bins; }

  void SetBins(std::vector<std::pair<double, double>> b) { Bins = b; }

  std::size_t numBins() { return Bins.size(); }

  std::vector<double> GetSystematics() { return sys; }

  std::vector<double> GetSystematicsError() { return sysError; }

  void SetSystematics(std::vector<double> b,
                      std::vector<double> bError = std::vector<double>());

  void SetSystematicsError(std::vector<double> bError);

  /// Anti particle systematics/
  std::vector<double> GetAntiSystematics() { return antiSys; }

  std::vector<double> GetAntiSystematicsError() { return antiSysError; }

  void SetAntiSystematics(std::vector<double> b,
                          std::vector<double> bError = std::vector<double>());

  void SetAntiSystematicsError(std::vector<double> bError);

  /** Count total systematics internally
   * Add systematics of track to total systematic error
   * @param charge
   * @param momentum
   */
  void AddToTotalError(int charge, double momentum);

  /// Total systematic uncertainty.
  /// The weighted mean of the uncertainties of the single tracks is calculated
  double GetTotalSystematics(bool useArithmeticMean = 0);

  /// Error of total systematic uncertainty.
  /// The weighted error of the uncertainties of the single tracks is
  /// calculated.
  double GetTotalSystematicsError(bool useArithmeticMean = 0);

  void Print() const;

  std::string GetTitle() { return title; }

  void SetTitle(std::string t) { title = t; }

protected:
  std::string title;

  /// invert e_mc/e_data-1 to e_data/e_mc-1
  static double inverse(double x);

  /// Calculate error for inversion e_mc/e_data-1 to e_data/e_mc-1
  static double inverseError(double x, double xErr);

  int findBin(double momentum) const;

  /// check if there is a bin defined that overlaps with [min,max]
  int findBin(double min, double max) const;

  void addBin(double min, double max);

  /// check for consistency
  bool check() const;

  std::vector<std::pair<double, double>> Bins;

  /// Data/MC difference in momentum bins for particle
  std::vector<double> sys;

  std::vector<double> sysError;

  /// Data/MC difference in momentum bins for anti-particle
  std::vector<double> antiSys;

  std::vector<double> antiSysError;

  std::vector<double> totalSys, totalSysError;
};

} // namespace Data
} // namespace ComPWA
#endif
