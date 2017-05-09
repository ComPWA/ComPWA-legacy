/*
 * MinuitResult.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: weidenka
 */

#include <numeric>
#include <cmath>
#include <cstdlib>

#include <boost/archive/xml_oarchive.hpp>

#include "Core/ProgressBar.hpp"
#include "Core/Logging.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

using namespace boost::log;

using ComPWA::Estimator::Estimator;

MinuitResult::MinuitResult()
    : calcInterference(0), initialLH(0), finalLH(0), trueLH(0) {}

MinuitResult::MinuitResult(std::shared_ptr<ControlParameter> esti,
                           ROOT::Minuit2::FunctionMinimum result)
    : calcInterference(0), initialLH(0), finalLH(0), trueLH(0) {
  est = std::static_pointer_cast<Estimator>(esti);
  _intens = est->GetIntensity();
  init(result);
}

void MinuitResult::setResult(std::shared_ptr<ControlParameter> esti,
                             ROOT::Minuit2::FunctionMinimum result) {
  est = std::static_pointer_cast<Estimator>(esti);
  _intens = est->GetIntensity();
  init(result);
}

void MinuitResult::init(ROOT::Minuit2::FunctionMinimum min) {
  ROOT::Minuit2::MnUserParameterState minState = min.UserState();

  if (minState.HasCovariance()) {
    ROOT::Minuit2::MnUserCovariance minuitCovMatrix = minState.Covariance();
    /* Size of Minuit covariance vector is given by dim*(dim+1)/2.
     * dim is the dimension of the covariance matrix.
     * The dimension can therefore be calculated as
     * dim = -0.5+-0.5 sqrt(8*size+1)
     */
    nFreeParameter = minuitCovMatrix.Nrow();
    globalCC = minState.GlobalCC().GlobalCC();
    cov = std::vector<std::vector<double>>(nFreeParameter,
                                           std::vector<double>(nFreeParameter));
    corr = std::vector<std::vector<double>>(
        nFreeParameter, std::vector<double>(nFreeParameter));
    for (unsigned i = 0; i < nFreeParameter; ++i)
      for (unsigned j = i; j < nFreeParameter; ++j) {
        cov.at(i).at(j) = minuitCovMatrix(j, i);
        cov.at(j).at(i) = minuitCovMatrix(j, i); // fill lower half
      }
    for (unsigned i = 0; i < nFreeParameter; ++i)
      for (unsigned j = i; j < nFreeParameter; ++j) {
        corr.at(i).at(j) =
            cov.at(i).at(j) / sqrt(cov.at(i).at(i) * cov.at(j).at(j));
        corr.at(j).at(i) = corr.at(i).at(j); // fill lower half
      }

  } else
    LOG(error) << "MinuitResult: no valid correlation matrix available!";
  initialLH = -1;
  finalLH = minState.Fval();
  edm = minState.Edm();
  isValid = min.IsValid();
  covPosDef = min.HasPosDefCovar();
  hasValidParameters = min.HasValidParameters();
  hasValidCov = min.HasValidCovariance();
  hasAccCov = min.HasAccurateCovar();
  hasReachedCallLimit = min.HasReachedCallLimit();
  edmAboveMax = min.IsAboveMaxEdm();
  hesseFailed = min.HesseFailed();
  errorDef = min.Up();
  nFcn = min.NFcn();

  return;
}

void MinuitResult::genSimpleOutput(std::ostream &out) {
  for (unsigned int o = 0; o < finalParameters.GetNDouble(); o++) {
    std::shared_ptr<DoubleParameter> outPar =
        finalParameters.GetDoubleParameter(o);
    out << outPar->GetValue() << " ";
    if (outPar->HasError())
      out << outPar->GetError() << " ";
  }
  out << "\n";

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
  out << "CPU Time : " << time / 60 << "min" << std::endl;
  out << std::setprecision(5) << std::endl;

  if (!hasValidParameters)
    out << "		*** NO VALID SET OF PARAMETERS! ***" << std::endl;
  if (printParam) {
    out << "PARAMETERS:" << std::endl;
    TableFormater *tableResult = new TableFormater(&out);
    PrintFitParameters(tableResult);
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
      PrintCovarianceMatrix(tableCov);
    }
    if (printCorrMatrix) {
      out << "CORRELATION MATRIX:" << std::endl;
      TableFormater *tableCorr = new TableFormater(&out);
      PrintCorrelationMatrix(tableCorr);
    }
  }
  out << "FIT FRACTIONS:" << std::endl;
  TableFormater tab(&out);
  PrintFitFractions(&tab);

  out << std::setprecision(10);
  out << "FinalLH: " << finalLH << std::endl;

  out << std::setprecision(5); // reset cout precision
  return;
}

void MinuitResult::createInterferenceTable(std::ostream &out,
                                           std::shared_ptr<AmpIntensity> amp) {
  out << "INTERFERENCE terms for " << amp->GetName() << ": " << std::endl;
  TableFormater *tableInterf = new TableFormater(&out);
  tableInterf->addColumn("Name 1", 15);
  tableInterf->addColumn("Name 2", 15);
  tableInterf->addColumn("Value", 15);
  tableInterf->header();
  double sumInfTerms = 0;
//  auto it = amp->GetResonanceItrFirst();
//  for (; it != amp->GetResonanceItrLast(); ++it) {
//    auto it2 = it;
//    for (; it2 != amp->GetResonanceItrLast(); ++it2) {
//      *tableInterf << (*it)->GetName();
//      *tableInterf << (*it2)->GetName();
//      double inf = amp->GetIntegralInterference(it, it2);
//      *tableInterf << inf;
//      sumInfTerms += inf;
//    }
//  }
  tableInterf->delim();
  *tableInterf << " "
               << "Sum: " << sumInfTerms;
  tableInterf->footer();
  out << std::endl;
}

void MinuitResult::PrintCorrelationMatrix(TableFormater *tableCorr) {
  if (!hasValidCov)
    return;
  tableCorr->addColumn(" ", 15);        // add empty first column
  tableCorr->addColumn("GlobalCC", 10); // global correlation coefficient

  // add columns in correlation matrix
  for (unsigned int o = 0; o < finalParameters.GetNDouble(); o++) {
    std::shared_ptr<DoubleParameter> ppp =
        finalParameters.GetDoubleParameter(o);
    if (ppp->IsFixed())
      continue;
    tableCorr->addColumn(ppp->GetName(), 15);
  }

  unsigned int n = 0;
  tableCorr->header();
  for (unsigned int o = 0; o < finalParameters.GetNDouble(); o++) {
    std::shared_ptr<DoubleParameter> ppp =
        finalParameters.GetDoubleParameter(o);
    if (ppp->IsFixed())
      continue;
    *tableCorr << ppp->GetName();
    *tableCorr << globalCC.at(n);
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

void MinuitResult::PrintCovarianceMatrix(TableFormater *tableCov) {
  if (!hasValidCov)
    return;
  tableCov->addColumn(" ", 17); // add empty first column
  // add columns first
  for (unsigned int o = 0; o < finalParameters.GetNDouble(); o++) {
    if (!finalParameters.GetDoubleParameter(o)->IsFixed())
      tableCov->addColumn(finalParameters.GetDoubleParameter(o)->GetName(), 17);
  }

  unsigned int n = 0;
  tableCov->header();
  for (unsigned int o = 0; o < finalParameters.GetNDouble(); o++) {
    std::shared_ptr<DoubleParameter> ppp =
        finalParameters.GetDoubleParameter(o);
    if (ppp->IsFixed())
      continue;
    *tableCov << ppp->GetName();
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

void MinuitResult::WriteXML(std::string filename) {
  std::ofstream ofs(filename);
  boost::archive::xml_oarchive oa(ofs);
  // TODO: Compile error due to serialization
  oa << boost::serialization::make_nvp("FitParameters", finalParameters);
  oa << boost::serialization::make_nvp("FitFractions", _fitFractions);
  ofs.close();
  return;
}

void MinuitResult::WriteTeX(std::string filename) {
  std::ofstream out(filename);
  TableFormater *tableResult = new TexTableFormater(&out);
  PrintFitParameters(tableResult);
  if (hasValidCov) {
    TableFormater *tableCov = new TexTableFormater(&out);
    PrintCovarianceMatrix(tableCov);
    TableFormater *tableCorr = new TexTableFormater(&out);
    PrintCorrelationMatrix(tableCorr);
  }
  TableFormater *fracTable = new TexTableFormater(&out);
  // calculate and print fractions if amplitude is set
  PrintFitFractions(fracTable);
  out.close();
  return;
}

bool MinuitResult::HasFailed() {
  bool failed = 0;
  if (!isValid)
    failed = 1;
  //	if(!covPosDef) failed=1;
  //	if(!hasValidParameters) failed=1;
  //	if(!hasValidCov) failed=1;
  //	if(!hasAccCov) failed=1;
  //	if(hasReachedCallLimit) failed=1;
  //	if(hesseFailed) failed=1;

  return failed;
}

} /* namespace Minuit2 */
} /* namespace Optimizer */
} /* namespace ComPWA */
