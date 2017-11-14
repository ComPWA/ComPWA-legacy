 
// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/FitResult.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {

void FitResult::WriteText(std::string filename) {
  std::ofstream myfile;
  myfile.open(filename, std::ios::app);
  genOutput(myfile);
  myfile.close();
  return;
};

void FitResult::WriteSimpleText(std::string filename) {
  std::ofstream myfile;
  myfile.open(filename);
  genSimpleOutput(myfile);
  myfile.close();
  return;
};

double FitResult::shiftAngle(double v) {
  double originalVal = v;
  double val = originalVal;
  while (val > M_PI)
    val -= 2 * M_PI;
  while (val < (-1) * M_PI)
    val += 2 * M_PI;
  if (val != originalVal)
    LOG(info) << "shiftAngle(): shifting parameter from " << originalVal
              << " to " << val << "!";
  return val;
}

void FitResult::genSimpleOutput(std::ostream &out) {
  for (unsigned int o = 0; o < finalParameters.GetNDouble(); o++) {
    std::shared_ptr<DoubleParameter> outPar =
        finalParameters.GetDoubleParameter(o);
    out << outPar->GetValue() << " " << outPar->GetError() << " ";
  }
  out << "\n";

  return;
}

void FitResult::SetFinalParameters(ParameterList &finPars) {
  finalParameters.DeepCopy(finPars);
}

void FitResult::SetTrueParameters(ParameterList &truePars) {
  trueParameters.DeepCopy(truePars);
}

void FitResult::SetFitFractions(ParameterList &list) {
  _fitFractions.DeepCopy(list);
}

void FitResult::Print(std::string opt) {
  std::stringstream s;
  genOutput(s, opt);
  std::string str = s.str();
  LOG(info) << str;
}

void FitResult::PrintFitParameters(TableFormater *tableResult) {
  bool printTrue = 0, printInitial = 0;
  if (trueParameters.GetNParameter())
    printTrue = 1;
  if (initialParameters.GetNParameter())
    printInitial = 1;

  // Column width for parameter with symmetric error
  unsigned int parErrorWidth = 22;

  // Do we have a parameter with assymetric errors?
  for (unsigned int o = 0; o < finalParameters.GetNDouble(); o++)
    if (finalParameters.GetDoubleParameter(o)->GetErrorType() ==
        ErrorType::ASYM)
      parErrorWidth = 33;

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

  for (unsigned int o = 0; o < finalParameters.GetNDouble(); o++) {
    std::shared_ptr<DoubleParameter> iniPar, outPar, truePar;
    try {
      outPar = finalParameters.GetDoubleParameter(o);
    } catch (BadParameter &bad) {
      LOG(error) << "FitResult::printFitParameters() | "
                    "can't access parameter of final parameter list!";
      throw;
    }
    if (printInitial) {
      try {
        iniPar = initialParameters.GetDoubleParameter(outPar->GetName());
      } catch (BadParameter &bad) {
        iniPar.reset();
      }
    }
    if (printTrue) {
      try {
        truePar = trueParameters.GetDoubleParameter(outPar->GetName());
      } catch (BadParameter &bad) {
        truePar.reset();
      }
    }

    ErrorType errorType = outPar->GetErrorType();
    bool isFixed = outPar->IsFixed();

    // Is parameter an angle?
    bool isAngle = 0;
    if (outPar->GetName().find("phase") != std::string::npos)
      isAngle = 1;
    // ... then shift the value to the domain [-pi;pi]
    if (isAngle && !isFixed) {
      outPar->SetValue(shiftAngle(outPar->GetValue()));
      if (printInitial)
        iniPar->SetValue(shiftAngle(iniPar->GetValue()));
      if (printTrue)
        truePar->SetValue(shiftAngle(truePar->GetValue()));
    }

    // Is parameter a magnitude?
    bool isMag = 0;
    if (outPar->GetName().find("mag") != std::string::npos)
      isMag = 1;
    // ... then make sure that it is positive
    if (isMag && !isFixed) {
      outPar->SetValue(std::fabs(outPar->GetValue()));
      if (printInitial && iniPar)
        iniPar->SetValue(std::fabs(iniPar->GetValue()));
      if (printTrue && truePar)
        truePar->SetValue(std::fabs(truePar->GetValue()));
    }

    // Print parameter name
    *tableResult << o << outPar->GetName();

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
      *tableResult << *outPar; // final value
    else
      *tableResult << " ";

    // Print true values
    if (printTrue) {
      if (truePar) {
        *tableResult << *truePar;
        double pull = (truePar->GetValue() - outPar->GetValue());
        // Shift pull by 2*pi if that reduces the deviation
        if (isAngle && !isFixed) {
          while (pull < 0 && pull < -M_PI)
            pull += 2 * M_PI;
          while (pull > 0 && pull > M_PI)
            pull -= 2 * M_PI;
        }
        if (outPar->HasError()) {
          if (errorType == ErrorType::ASYM && pull < 0)
            pull /= outPar->GetErrorLow();
          else if (errorType == ErrorType::ASYM && pull > 0)
            pull /= outPar->GetErrorHigh();
          else
            pull /= outPar->GetError();
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

void FitResult::PrintFitFractions(TableFormater *fracTable) {
  LOG(info) << " FitResult::printFitFractions() | "
               "Calculating fit fractions!";

  double sum = 0, sumErrorSq = 0;

  fracTable->Reset();

  // print matrix
  fracTable->addColumn("Fit fractions [%]", 40); // add empty first column
  fracTable->addColumn("Fraction", 15);          // add empty first column
  fracTable->addColumn("Error", 15);             // add empty first column
  fracTable->header();
  for (unsigned int i = 0; i < _fitFractions.GetNDouble(); ++i) {
    std::shared_ptr<DoubleParameter> tmpPar =
        _fitFractions.GetDoubleParameter(i);
    std::string resName = tmpPar->GetName();

    *fracTable << resName << tmpPar->GetValue();
    try {
      *fracTable << tmpPar->GetError(); // assume symmetric errors here
    } catch (std::exception &ex) {
      *fracTable << 0.0;
    }

    sum += tmpPar->GetValue();
    sumErrorSq += tmpPar->GetError() * tmpPar->GetError();
  }
  fracTable->delim();
  *fracTable << "Total" << sum << sqrt(sumErrorSq);
  fracTable->footer();
  sumFractions = sum;
  sumFractionsError = sqrt(sumErrorSq);

  return;
}

} // namespace 
