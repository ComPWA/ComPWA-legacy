// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Optimizer/Geneva/GenevaResult.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Geneva {

void GenevaResult::init(std::shared_ptr<Gem::Geneva::GStartIndividual> min) {
  //	finalLH = min->getBestCase();
  finalLH =
      std::get<1>(min->getBestKnownPrimaryFitness()); // TODO: Check if correct
                                                      // way to get result
  min->getPar(FinalParameters);
  // ToDO: extract more info
  return;
}

void GenevaResult::genOutput(std::ostream &out, std::string opt) {
  bool printParam = 1, printCorrMatrix = 1, printCovMatrix = 1;
  if (opt == "P") { // print only parameters
    printCorrMatrix = 0;
    printCovMatrix = 0;
  }
  out << std::endl;
  out << "--------------MINUIT2 FIT RESULT----------------" << std::endl;
  if (!isValid)
    out << "		*** MINIMUM NOT VALID! ***" << std::endl;
  out << std::setprecision(10);
  out << "Initial Likelihood: " << initialLH << std::endl;
  out << "Final Likelihood: " << finalLH << std::endl;
  if (trueLH)
    out << "True Likelihood: " << trueLH << std::endl;

  out << "Estimated distance to minimumn: " << edm << std::endl;
  if (edmAboveMax)
    out << "		*** EDM IS ABOVE MAXIMUM! ***" << std::endl;
  out << "Error definition: " << errorDef << std::endl;
  out << "Number of calls: " << nFcn << std::endl;
  if (hasReachedCallLimit)
    out << "		*** LIMIT OF MAX CALLS REACHED! ***" << std::endl;
  out << "CPU Time : " << time() / 60 << "min" << std::endl;
  out << std::setprecision(5) << std::endl;

  if (!hasValidParameters)
    out << "		*** NO VALID SET OF PARAMETERS! ***" << std::endl;
  if (printParam) {
    out << "PARAMETERS:" << std::endl;
    TableFormater *tableResult = new TableFormater(&out);
    printFitParameters(tableResult);
  }

  if (!hasValidCov)
    out << "		*** COVARIANCE MATRIX NOT VALID! ***" << std::endl;
  if (!hasAccCov)
    out << "		*** COVARIANCE MATRIX NOT ACCURATE! ***" << std::endl;
  if (!covPosDef)
    out << "		*** COVARIANCE MATRIX NOT POSITIVE DEFINITE! ***"
        << std::endl;
  if (hesseFailed)
    out << "		*** HESSE FAILED! ***" << std::endl;
  if (hasValidCov) {
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
  out << "FinalLH: " << finalLH << std::endl;

  out << std::setprecision(5); // reset cout precision
  return;
}

void GenevaResult::printCorrelationMatrix(TableFormater *tableCorr) {
  if (!hasValidCov)
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
      *tableCorr << globalCC.at(n);
    } catch (...) {
      *tableCorr << "?";
    }

    for (unsigned int t = 0; t < corr.size(); t++) {
      if (n >= corr.at(0).size()) {
        *tableCorr << " ";
        continue;
      }
      if (t >= n)
        *tableCorr << corr.at(n).at(t);
      else
        *tableCorr << "";
    }
    n++;
  }
  tableCorr->footer();
  return;
}

void GenevaResult::printCovarianceMatrix(TableFormater *tableCov) {
  if (!hasValidCov)
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
    for (unsigned int t = 0; t < cov.size(); t++) {
      if (n >= cov.at(0).size()) {
        *tableCov << " ";
        continue;
      }
      if (t >= n)
        *tableCov << cov.at(n).at(t);
      else
        *tableCov << "";
    }
    n++;
  }
  tableCov->footer();
  return;
}

} // namespace Geneva
} // namespace Optimizer
} // namespace ComPWA
