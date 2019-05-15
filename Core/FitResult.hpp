// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_FITRESULT_HPP_
#define COMPWA_FITRESULT_HPP_

#include <unordered_map>

#include "Core/ParameterList.hpp"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/serialization.hpp>

namespace ComPWA {

struct FitResult {
  // these properties are separate members, since they should always be set and
  // are used most commonly
  ParameterList InitialParameters;
  ParameterList FinalParameters;

  double InitialEstimatorValue = 0.0;
  double FinalEstimatorValue = 0.0;
  double ElapsedTimeInSeconds = 0.0;

  /// optional properties
  std::unordered_map<std::string, double> DoubleProperties;
  std::unordered_map<std::string, int> IntProperties;
  std::unordered_map<std::string, ParameterList> ParameterListProperties;
  std::unordered_map<std::string, std::vector<double>> DoubleListProperties;
  std::unordered_map<std::string, std::vector<std::vector<double>>>
      MatrixProperties;

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    ar &BOOST_SERIALIZATION_NVP(ElapsedTimeInSeconds);
    ar &BOOST_SERIALIZATION_NVP(InitialEstimatorValue);
    ar &BOOST_SERIALIZATION_NVP(FinalEstimatorValue);
    ar &BOOST_SERIALIZATION_NVP(InitialParameters);
    ar &BOOST_SERIALIZATION_NVP(FinalParameters);
    //ar << DoubleProperties;
    //ar &BOOST_SERIALIZATION_NVP(DoubleProperties);
    //ar &BOOST_SERIALIZATION_NVP(IntProperties);
    //ar &BOOST_SERIALIZATION_NVP(MatrixProperties);
  }
};

void printFitResult(const FitResult &result) {

}

/*
void printFitResult(std::ostream &out, std::string opt) {
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



void printFitResult(const FitResult &fitresult) {
	//TableFormater *tableResult) {

	bool printTrue = 0, printInitial = 0;
  if (TrueParameters.numParameters())
    printTrue = 1;
  if (InitialParameters.numParameters())
    printInitial = 1;

  // Column width for parameter with symmetric error
  size_t parErrorWidth = 22;

  // Do we have a parameter with assymetric errors?
  //  for (unsigned int o = 0; o < finalParameters.GetNDouble(); o++)
  for (auto p : FinalParameters.doubleParameters()) {
    if (p->errorType() == ErrorType::ASYM)
      parErrorWidth = 33;
  }

  tableResult->addColumn("Nr");
  tableResult->addColumn("Name", 30);
  if (printInitial)
    tableResult->addColumn("Initial Value", parErrorWidth);
  tableResult->addColumn("Final Value", parErrorWidth);
  if (printTrue)
    tableResult->addColumn("True Value", 10);
  if (printTrue)
    tableResult->addColumn("Deviation", 9);
  tableResult->header();

  size_t parameterId = 0;
  for (auto p : FinalParameters.doubleParameters()) {
    //    for (unsigned int o = 0; o < finalParameters.GetNDouble(); o++) {
    std::shared_ptr<FitParameter> iniPar, truePar;
    std::string name = p->name();

    if (printInitial) {
      try {
        iniPar = FindParameter(p->name(), InitialParameters);
      } catch (BadParameter &bad) {
        iniPar.reset();
      }
    }
    if (printTrue) {
      try {
        truePar = FindParameter(p->name(), TrueParameters);
      } catch (BadParameter &bad) {
        truePar.reset();
      }
    }

    ErrorType errorType = p->errorType();
    bool isFixed = p->isFixed();

    // Is parameter an angle?
    bool isAngle = 0;
    if (p->name().find("phase") != std::string::npos)
      isAngle = 1;
    // ... then shift the value to the domain [-pi;pi]
    if (isAngle && !isFixed) {
      p->setValue(shiftAngle(p->value()));
      if (printInitial)
        iniPar->setValue(shiftAngle(iniPar->value()));
      if (printTrue)
        truePar->setValue(shiftAngle(truePar->value()));
    }

    // Is parameter a magnitude?
    bool isMag = 0;
    if (p->name().find("mag") != std::string::npos)
      isMag = 1;
    // ... then make sure that it is positive
    if (isMag && !isFixed) {
      p->setValue(std::abs(p->value()));
      if (printInitial && iniPar)
        iniPar->setValue(std::abs(iniPar->value()));
      if (printTrue && truePar)
        truePar->setValue(std::abs(truePar->value()));
    }

    // Print parameter name
    *tableResult << parameterId << p->name();
    parameterId++;

    // Print initial values
    if (printInitial) {
      if (iniPar) {
        (*tableResult) << *iniPar; // |nr.| name| inital value|
      } else {
        (*tableResult) << " ";
      }
    }

    // Print final value
    if (!isFixed)
      *tableResult << *p; // final value
    else
      *tableResult << " ";

    // Print true values
    if (printTrue) {
      if (truePar) {
        *tableResult << *truePar;
        double pull = (truePar->value() - p->value());
        // Shift pull by 2*pi if that reduces the deviation
        if (isAngle && !isFixed) {
          while (pull < 0 && pull < -M_PI)
            pull += 2 * M_PI;
          while (pull > 0 && pull > M_PI)
            pull -= 2 * M_PI;
        }
        if (p->hasError()) {
          if (errorType == ErrorType::ASYM && pull < 0)
            pull /= p->error().first;
          else if (errorType == ErrorType::ASYM && pull > 0)
            pull /= p->error().second;
          else
            pull /= p->avgError();
        }
        if (!std::isnan(pull))
          (*tableResult) << pull;
        else
          (*tableResult) << " ";
      } else {
        (*tableResult) << " "
                       << " ";
      }
    }
  }
  tableResult->footer();
  return;
}

void printFitFractions(TableFormater *fracTable) {
  LOG(INFO) << " FitResult::printFitFractions() | "
               "Calculating fit fractions!";

  double sum = 0, sumErrorSq = 0;

  fracTable->reset();

  // print matrix
  fracTable->addColumn("Fit fractions [%]", 40); // add empty first column
  fracTable->addColumn("Fraction", 15);          // add empty first column
  fracTable->addColumn("Error", 15);             // add empty first column
  fracTable->header();
  for( auto f :FitFractions.doubleParameters() ){
    std::string resName = f->name();
    *fracTable << resName << f->value();
    try {
      *fracTable << f->avgError(); // assume symmetric errors here
    } catch (std::exception &ex) {
      *fracTable << 0.0;
    }
    sum += f->value();
    sumErrorSq += f->error().first * f->error().first;
  }

  fracTable->delim();
  *fracTable << "Total" << sum << sqrt(sumErrorSq);
  fracTable->footer();
  SumFractions = sum;
  SumFractionsError = sqrt(sumErrorSq);

  return;
}


void createInterferenceTable(std::ostream &out) {
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

void printCorrelationMatrix(TableFormater *tableCorr) {
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

void printCovarianceMatrix(TableFormater *tableCov) {
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

void writeXML(std::string filename) {
  std::ofstream ofs(filename);
  boost::archive::xml_oarchive oa(ofs);
  // TODO: Compile error due to serialization
  oa << boost::serialization::make_nvp("FitParameters", FinalParameters);
  oa << boost::serialization::make_nvp("FitFractions", FitFractions);
  ofs.close();
  return;
}*/

} // namespace ComPWA

#endif
