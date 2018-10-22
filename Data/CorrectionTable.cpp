// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.
#include "Data/CorrectionTable.hpp"

#include <ostream>
#include <string>

#include "Core/Logging.hpp"
#include "Core/TableFormater.hpp"

namespace ComPWA {
namespace Data {

void CorrectionTable::Print() const {
  if (sys.size() == 0 && antiSys.size() == 0)
    return; // don't print if empty
  std::stringstream out;
  ComPWA::TableFormater table(&out);
  table.addColumn("Bin", 12);             // add empty first column
  table.addColumn("Particle", 20);        // global correlation coefficient
  table.addColumn("anti-Particle", 20);   // global correlation coefficient
  table.addColumn("Neutral/Average", 20); // global correlation coefficient
  table.header();
  for (unsigned int i = 0; i < Bins.size(); i++) {
    std::stringstream strBin;
    strBin << Bins.at(i).first << " - " << Bins.at(i).second;
    table << strBin.str();
    std::stringstream strCor;
    strCor << GetValue(+1, Bins.at(i).second);
    //				<< " +- " <<GetError(+1,Bins.at(i).second);
    table << strCor.str();
    std::stringstream strAntiCor;
    strAntiCor << GetValue(-1, Bins.at(i).second);
    //				<< " +- " <<GetError(-1,Bins.at(i).second);
    table << strAntiCor.str();
    std::stringstream strAvg;
    strAvg << GetValue(0, Bins.at(i).second);
    //				<< " +- " <<GetError(0,Bins.at(i).second);
    table << strAvg.str();
  }
  table.footer();
  LOG(INFO) << "CorrectionTable::Print() | " << title << std::endl << out.str();
}
void CorrectionTable::SetSystematics(std::vector<double> b,
                                     std::vector<double> bError) {
  sys = b;
  if (bError.size() == b.size()) // is error vector valid?
    sysError = bError;
  else // otherwise set errors to 0
    sysError = std::vector<double>(Bins.size(), 0.);
  if (!antiSys.size()) { // if anti-particle systematics is empty set is to the
                         // same values
    antiSys = b;
    antiSysError = bError;
  }
  return;
}
void CorrectionTable::SetAntiSystematics(std::vector<double> b,
                                         std::vector<double> bError) {
  antiSys = b;
  if (bError.size() == b.size()) // is error vector valid?
    antiSysError = bError;
  else // otherwise set errors to 0
    antiSysError = std::vector<double>(Bins.size(), 0.);
  return;
}
void CorrectionTable::SetSystematicsError(std::vector<double> bError) {
  sysError = bError;
  return;
}
void CorrectionTable::SetAntiSystematicsError(std::vector<double> bError) {
  antiSysError = bError;
  return;
}
void CorrectionTable::AddInv(double binMin, double binMax, double s,
                             double sError, double antiS, double antiSerror) {
  Add(binMin, binMax, inverse(s), inverseError(s, sError), inverse(antiS),
      inverseError(antiS, antiSerror));
}
void CorrectionTable::Add(double binMin, double binMax, double s, double sError,
                          double antiS, double antiSerror) {
  addBin(binMin, binMax);
  sys.push_back(s);
  sysError.push_back(sError);
  if (antiS == -999) { // no efficiency is given for anti-particles
    antiSys.push_back(s);
    antiSysError.push_back(sError);
  } else {
    antiSys.push_back(antiS);
    antiSysError.push_back(antiSerror);
  }
}
double CorrectionTable::GetInvValue(int charge, double momentum) const {
  return inverse(GetValue(charge, momentum));
}
double CorrectionTable::GetInvError(int charge, double momentum) const {
  return inverseError(GetValue(charge, momentum), GetError(charge, momentum));
}
double CorrectionTable::GetValue(int charge, double momentum) const {
  int binNumber = findBin(momentum);
  if (binNumber < 0)
    throw std::runtime_error(
        "momentumSys::value() no bin found for momentum value " +
        std::to_string((long double)momentum) + "!");

  double val = 0.0;
  if (charge > 0)
    val = sys.at(binNumber);
  else if (charge < 0)
    val = antiSys.at(binNumber);
  else if (charge == 0) { // charge unknown -> average value
    val = (antiSys.at(binNumber) + sys.at(binNumber)) / 2;
  }
  return val;
}
double CorrectionTable::GetError(int charge, double momentum) const {
  int binNumber = findBin(momentum);
  if (binNumber < 0)
    throw std::runtime_error(
        "momentumSys::GetError() no bin found for momentum value " +
        std::to_string((long double)momentum) + "!");

  double err = 0.0;
  if (charge > 0) // D0->K0bar K+K-
    err = sysError.at(binNumber);
  else if (charge < 0) // D0->K0bar K+K-
    err = antiSysError.at(binNumber);
  else if (charge == 0) { // charge unknown -> average value
    err = (0.5 *
           std::sqrt(antiSysError.at(binNumber) * antiSysError.at(binNumber) +
                     sysError.at(binNumber) * sysError.at(binNumber)));
  }
  return err;
}

void CorrectionTable::AddToTotalError(int charge, double momentum) {
  try {
    totalSys.push_back(GetValue(charge, momentum));
    totalSysError.push_back(GetError(charge, momentum));
  } catch (const std::exception &e) {
    LOG(INFO) << "momentumSys::AddToTotalError()  can't add "
                 "systematics for momentum "
              << momentum << " to total systematics. We skip that track!";
  }
}

double CorrectionTable::GetTotalSystematics(bool useArithmeticMean) {
  double mean = 0;
  double sumSigma = 0;
  // use arithmetic mean
  if (useArithmeticMean) {
    for (unsigned int i = 0; i < totalSys.size(); i++)
      mean += totalSys.at(i);
    return mean / totalSys.size();
  } else {
    // use weighted mean
    for (unsigned int i = 0; i < totalSys.size(); i++) {
      double tmp = 1 / (totalSysError.at(i) * totalSysError.at(i));
      mean += tmp * totalSys.at(i);
      sumSigma += tmp;
    }
    mean /= sumSigma;
    return mean;
  }
}

double CorrectionTable::GetTotalSystematicsError(bool useArithmeticMean) {
  double sumSigma = 0; // inverse sum over uncertainties
  double mean = GetTotalSystematics(useArithmeticMean);
  if (useArithmeticMean) {
    for (unsigned int i = 0; i < totalSys.size(); i++)
      sumSigma += (totalSys.at(i) - mean) * (totalSys.at(i) - mean);
    return sumSigma = sqrt(sumSigma / totalSys.size());
  } else {
    for (unsigned int i = 0; i < totalSysError.size(); i++)
      sumSigma += 1 / (totalSysError.at(i) * totalSysError.at(i));
    return sqrt(1 / sumSigma);
  }
}

double CorrectionTable::inverse(double x) {
  if (x == -999 || x == -1)
    return -999;
  return 1 / (x + 1) - 1;
}

double CorrectionTable::inverseError(double x, double xErr) {
  if (x == -999 || x == -1)
    return -999;
  return std::fabs(-1 / ((x + 1) * (x + 1)) * xErr);
}

int CorrectionTable::findBin(double momentum) const {
  for (unsigned int i = 0; i < Bins.size(); i++) {
    if (momentum <= Bins.at(i).second && momentum > Bins.at(i).first) {
      return i;
    }
  }
  return -1; // no bin found
}

int CorrectionTable::findBin(double min, double max) const {
  for (unsigned int i = 0; i < Bins.size(); i++) {
    if ((min < Bins.at(i).second && min >= Bins.at(i).first) ||
        (max <= Bins.at(i).second && max > Bins.at(i).first))
      return i;
  }
  return -1; // bin doesn't exist
}

void CorrectionTable::addBin(double min, double max) {
  if (findBin(min, max) < 0) // is a bin already defined?
    Bins.push_back(std::make_pair(min, max));
  else
    throw std::runtime_error("CorrectionTable::addBin() in range [" +
                             std::to_string((long double)min) + "," +
                             std::to_string((long double)max) +
                             "] a bin is already defined!");
}

bool CorrectionTable::check() const {
  unsigned int s = Bins.size();
  if (sys.size() != s || sysError.size() != s || antiSys.size() != s ||
      antiSysError.size() != s)
    return 0;
  return 1;
}

} // namespace Data
} // namespace ComPWA
