// Copyright (c) 2014, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
#include <cstdlib>
#include <numeric>

#include <Minuit2/MnUserParameterState.h>
#include <boost/archive/xml_oarchive.hpp>

#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"

using namespace ComPWA::Optimizer::Minuit2;

MinuitResult::MinuitResult()
    : est(nullptr), CalcInterference(0), InitialLH(0), FinalLH(0), TrueLH(0) {}

MinuitResult::MinuitResult(std::shared_ptr<ComPWA::Estimator::Estimator> esti,
                           ROOT::Minuit2::FunctionMinimum result)
    : est(esti), CalcInterference(0), InitialLH(0), FinalLH(0), TrueLH(0) {
  init(result);
}

void MinuitResult::setResult(std::shared_ptr<ComPWA::Estimator::Estimator> esti,
                             ROOT::Minuit2::FunctionMinimum result) {
  est = esti;
  init(result);
}

void MinuitResult::init(ROOT::Minuit2::FunctionMinimum min) {

  ROOT::Minuit2::MnUserParameterState minState = min.UserState();
  NumFreeParameter = minState.Parameters().Trafo().VariableParameters();

  if (minState.HasCovariance()) {
    ROOT::Minuit2::MnUserCovariance minuitCovMatrix = minState.Covariance();
    // Size of Minuit covariance vector is given by dim*(dim+1)/2.
    // dim is the dimension of the covariance matrix.
    // The dimension can therefore be calculated as
    // dim = -0.5+-0.5 sqrt(8*size+1)
    assert(minuitCovMatrix.Nrow() == NumFreeParameter);
    Cov = std::vector<std::vector<double>>(
        NumFreeParameter, std::vector<double>(NumFreeParameter));
    Corr = std::vector<std::vector<double>>(
        NumFreeParameter, std::vector<double>(NumFreeParameter));
    for (unsigned i = 0; i < NumFreeParameter; ++i)
      for (unsigned j = i; j < NumFreeParameter; ++j) {
        Cov.at(i).at(j) = minuitCovMatrix(j, i);
        Cov.at(j).at(i) = minuitCovMatrix(j, i); // fill lower half
      }
    for (unsigned i = 0; i < NumFreeParameter; ++i)
      for (unsigned j = i; j < NumFreeParameter; ++j) {
        Corr.at(i).at(j) =
            Cov.at(i).at(j) / sqrt(Cov.at(i).at(i) * Cov.at(j).at(j));
        Corr.at(j).at(i) = Corr.at(i).at(j); // fill lower half
      }

  } else {
    LOG(ERROR) << "MinuitResult: no valid correlation matrix available!";
    // Initialize empty covariance matrix with correct dimensions
    Cov = std::vector<std::vector<double>>(
        NumFreeParameter, std::vector<double>(NumFreeParameter, 0.));
    // Initialize with default values(?)
    //    for(unsigned int t = 0; t < n; t++)
    //      Cov[t][t] = 1;
  }
  if (minState.HasGlobalCC()) {
    GlobalCC = minState.GlobalCC().GlobalCC();
  } else {
    GlobalCC = std::vector<double>(NumFreeParameter, 0);
    LOG(ERROR) << "MinuitResult: no valid global correlation available!";
  }
  InitialLH = -1;
  FinalLH = minState.Fval();
  Edm = minState.Edm();
  IsValid = min.IsValid();
  CovPosDef = min.HasPosDefCovar();
  HasValidParameters = min.HasValidParameters();
  HasValidCov = min.HasValidCovariance();
  HasAccCov = min.HasAccurateCovar();
  HasReachedCallLimit = min.HasReachedCallLimit();
  EdmAboveMax = min.IsAboveMaxEdm();
  HesseFailed = min.HesseFailed();
  ErrorDef = min.Up();
  NFcn = min.NFcn();

  return;
}

void MinuitResult::genOutput(std::ostream &out, std::string opt) {
  bool printParam = 1, printCorrMatrix = 1, printCovMatrix = 1;
  if (opt == "P") { // print only parameters
    printCorrMatrix = 0;
    printCovMatrix = 0;
  }
  out << std::endl;
  out << "--------------MINUIT2 FIT RESULT----------------" << std::endl;
  if (!IsValid)
    out << "    *** MINIMUM NOT VALID! ***" << std::endl;
  out << std::setprecision(10);
  out << "Initial Likelihood: " << InitialLH << std::endl;
  out << "Final Likelihood: " << FinalLH << std::endl;
  if (TrueLH)
    out << "True Likelihood: " << TrueLH << std::endl;

  out << "Estimated distance to minimumn: " << Edm << std::endl;
  if (EdmAboveMax)
    out << "    *** EDM IS ABOVE MAXIMUM! ***" << std::endl;
  out << "Error definition: " << ErrorDef << std::endl;
  out << "Number of calls: " << NFcn << std::endl;
  if (HasReachedCallLimit)
    out << "    *** LIMIT OF MAX CALLS REACHED! ***" << std::endl;
  out << "CPU Time : " << Time / 60 << "min" << std::endl;
  out << std::setprecision(5) << std::endl;

  if (!HasValidParameters)
    out << "    *** NO VALID SET OF PARAMETERS! ***" << std::endl;
  if (printParam) {
    out << "PARAMETERS:" << std::endl;
    TableFormater *tableResult = new TableFormater(&out);
    printFitParameters(tableResult);
  }

  if (!HasValidCov)
    out << "    *** COVARIANCE MATRIX NOT VALID! ***" << std::endl;
  if (!HasAccCov)
    out << "    *** COVARIANCE MATRIX NOT ACCURATE! ***" << std::endl;
  if (!CovPosDef)
    out << "    *** COVARIANCE MATRIX NOT POSITIVE DEFINITE! ***" << std::endl;
  if (HesseFailed)
    out << "    *** HESSE FAILED! ***" << std::endl;
  if (HasValidCov) {
    if (printCovMatrix) {
      out << "COVARIANCE MATRIX:" << std::endl;
      TableFormater *tableCov = new TableFormater(&out);
      printCovarianceMatrix(tableCov);
    }
    if (printCorrMatrix) {
      out << "CORRELATION MATRIX:" << std::endl;
      TableFormater *tableCorr = new TableFormater(&out);
      printCorrelationMatrix(tableCorr);
    }
  }
  out << "FIT FRACTIONS:" << std::endl;
  TableFormater tab(&out);
  printFitFractions(&tab);

  out << std::setprecision(10);
  out << "FinalLH: " << FinalLH << std::endl;

  out << std::setprecision(5); // reset cout precision
  return;
}

void MinuitResult::createInterferenceTable(std::ostream &out,
                                           std::shared_ptr<Intensity> amp) {
  //out << "INTERFERENCE terms for " << amp->getName() << ": " << std::endl;
  TableFormater *tableInterf = new TableFormater(&out);
  tableInterf->addColumn("Name 1", 15);
  tableInterf->addColumn("Name 2", 15);
  tableInterf->addColumn("Value", 15);
  tableInterf->header();
  double sumInfTerms = 0;
  tableInterf->delim();
  *tableInterf << " "
               << "Sum: " << sumInfTerms;
  tableInterf->footer();
  out << std::endl;
}

void MinuitResult::printCorrelationMatrix(TableFormater *tableCorr) {
  if (!HasValidCov)
    return;
  tableCorr->addColumn(" ", 15);        // add empty first column
  tableCorr->addColumn("GlobalCC", 10); // global correlation coefficient

  // add columns in correlation matrix
  for (auto p : FinalParameters.doubleParameters()) {
    if (p->isFixed())
      continue;
    tableCorr->addColumn(p->name(), 15);
  }

  unsigned int n = 0;
  tableCorr->header();
  for (auto p : FinalParameters.doubleParameters()) {
    if (p->isFixed())
      continue;
    *tableCorr << p->name();
    try {
      *tableCorr << GlobalCC.at(n);
    } catch (...) {
      *tableCorr << "?";
    }

    for (unsigned int t = 0; t < Corr.size(); t++) {
      if (n >= Corr.at(0).size()) {
        *tableCorr << " ";
        continue;
      }
      if (t >= n)
        *tableCorr << Corr.at(n).at(t);
      else
        *tableCorr << "";
    }
    n++;
  }
  tableCorr->footer();
  return;
}

void MinuitResult::printCovarianceMatrix(TableFormater *tableCov) {
  if (!HasValidCov)
    return;
  // Create table structure first
  tableCov->addColumn(" ", 17); // add empty first column
                                // add columns first
  for (auto p : FinalParameters.doubleParameters()) {
    if (p->isFixed())
      continue;
    tableCov->addColumn(p->name(), 17);
  }

  // Fill table
  unsigned int n = 0;
  tableCov->header();
  for (auto p : FinalParameters.doubleParameters()) {
    if (p->isFixed())
      continue;

    *tableCov << p->name();
    for (unsigned int t = 0; t < Cov.size(); t++) {
      if (n >= Cov.at(0).size()) {
        *tableCov << " ";
        continue;
      }
      if (t >= n)
        *tableCov << Cov.at(n).at(t);
      else
        *tableCov << "";
    }
    n++;
  }
  tableCov->footer();
  return;
}

void MinuitResult::writeXML(std::string filename) {
  std::ofstream ofs(filename);
  boost::archive::xml_oarchive oa(ofs);
  // TODO: Compile error due to serialization
  oa << boost::serialization::make_nvp("FitParameters", FinalParameters);
  oa << boost::serialization::make_nvp("FitFractions", FitFractions);
  ofs.close();
  return;
}

void MinuitResult::writeTeX(std::string filename) {
  std::ofstream out(filename);
  TableFormater *tableResult = new TexTableFormater(&out);
  printFitParameters(tableResult);
  if (HasValidCov) {
    TableFormater *tableCov = new TexTableFormater(&out);
    printCovarianceMatrix(tableCov);
    TableFormater *tableCorr = new TexTableFormater(&out);
    printCorrelationMatrix(tableCorr);
  }
  TableFormater *fracTable = new TexTableFormater(&out);
  // calculate and print fractions if amplitude is set
  printFitFractions(fracTable);
  out.close();
  return;
}
