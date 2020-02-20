
// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/FitResult.hpp"

#include <iomanip>

#include "Core/Logging.hpp"
#include "Core/TableFormatter.hpp"
#include "Core/Utils.hpp"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
namespace ComPWA {

std::string makeFitParameterString(ComPWA::FitParameter<double> p) {
  std::stringstream ss;
  ss << std::setprecision(10);
  ss << p.Value;
  ss << std::setprecision(5);
  if (p.IsFixed) {
    ss << " (fixed)";
  } else if (p.Error.first == p.Error.second)
    ss << " +- " << p.Error.first;
  else {
    ss << " [" << p.Error.first << ", " << p.Error.second << "]";
  }
  return ss.str();
}

std::ostream &operator<<(std::ostream &os, const FitResult &Result) {
  auto OldPrecision(os.precision());

  using ComPWA::TableFormatter;

  os << "\n--------------FIT RESULT INFOS----------------\n";
  // ----------------- general info -----------------
  os << "initial estimator value: " << Result.InitialEstimatorValue << "\n";
  os << "final estimator value: " << Result.FinalEstimatorValue << "\n";
  os << "duration of fit (in seconds): " << Result.FitDuration.count() << "\n";

  // ----------------- fit parameters -----------------
  TableFormatter FitParametersFormatter(&os);
  // Column width for parameter with symmetric error
  size_t ParErrorWidth = 30;

  // Any parameter with asymmetric errors?
  for (auto x : Result.FinalParameters) {
    if (x.Error.first != x.Error.second) {
      ParErrorWidth = 40;
      break;
    }
  }

  FitParametersFormatter.addColumn("Nr");
  FitParametersFormatter.addColumn("Name", 30);
  FitParametersFormatter.addColumn("Initial Value", 15);
  FitParametersFormatter.addColumn("Final Value", ParErrorWidth);

  FitParametersFormatter.header();

  os << std::setprecision(10);
  size_t parameterId = 0;
  for (auto p : Result.FinalParameters) {
    auto res = std::find_if(Result.InitialParameters.begin(),
                            Result.InitialParameters.end(),
                            [&p](const ComPWA::FitParameter<double> &x) {
                              return p.Name == x.Name;
                            });

    // Is parameter an angle?
    bool IsAngle(false);
    if (p.Name.find("phase") != std::string::npos)
      IsAngle = true;

    // Print parameter name
    FitParametersFormatter << parameterId << p.Name;
    parameterId++;

    // Print initial values
    if (res != Result.InitialParameters.end()) {
      if (IsAngle)
        FitParametersFormatter << ComPWA::Utils::shiftAngle(res->Value);
      else
        FitParametersFormatter << res->Value;
    } else {
      LOG(WARNING)
          << "FitResult operator<<(): could not find initial parameter. "
             "FitResult corrupted, skipping initial value";
      FitParametersFormatter << " "; // print blank into that table cell
    }

    // Print final value
    if (IsAngle)
      p.Value = ComPWA::Utils::shiftAngle(p.Value);
    FitParametersFormatter << makeFitParameterString(p);
  }

  FitParametersFormatter.footer();

  os << std::setprecision(5);
  // ----------------- covariance matrix -----------------
  size_t NRows(Result.CovarianceMatrix.size());
  if (0 < NRows) {
    bool CovarianceValid(true);
    for (auto x : Result.CovarianceMatrix) {
      if (x.size() != NRows) {
        CovarianceValid = false;
        LOG(WARNING)
            << "FitResult operator<<(): Covariance is not a square matrix!";
        break;
      }
    }
    if (CovarianceValid) {
      TableFormatter CovarianceFormatter(&os);

      // Create table structure first
      CovarianceFormatter.addColumn(" ", 17); // add empty first column
                                              // add columns first
      for (auto p : Result.FinalParameters) {
        if (p.IsFixed)
          continue;
        CovarianceFormatter.addColumn(p.Name, 17);
      }

      // Fill table
      unsigned int n = 0;
      CovarianceFormatter.header();
      for (auto p : Result.FinalParameters) {
        if (p.IsFixed)
          continue;

        CovarianceFormatter << p.Name;
        for (auto val : Result.CovarianceMatrix.at(n)) {
          CovarianceFormatter << val;
        }
        n++;
      }
      CovarianceFormatter.footer();
    }
  }
  os << std::setprecision(OldPrecision); // reset os precision

  return os;
}

void FitResult::write(std::string filename) const {
  std::ofstream ofs(filename);
  boost::archive::xml_oarchive oa(ofs);
  oa << boost::serialization::make_nvp("FitResult", *this);
}

FitResult load(std::string filename) {
  FitResult Result;
  std::ifstream ifs(filename);
  assert(ifs.good());
  boost::archive::xml_iarchive ia(ifs);
  ia >> boost::serialization::make_nvp("FitResult", Result);
  return Result;
}

void initializeWithFitResult(ComPWA::Intensity &Intens,
                             ComPWA::FitResult Result) {
  auto Params = Intens.getParameters();
  std::vector<double> ParameterValues;
  for (auto p : Params) {
    auto found = std::find_if(Result.FinalParameters.begin(),
                              Result.FinalParameters.end(),
                              [&p](auto const &x) { return p.Name == x.Name; });
    if (Result.FinalParameters.end() == found) {
      LOG(ERROR) << "initializeWithFitResult(): Could not find parameter "
                 << p.Name
                 << " in fit result! Intensity is not fully initialized!"
                 << "\nPlease make sure the supplied fit result "
                    "is compatible with this Intensity.";
    } else {
      ParameterValues.push_back(found->Value);
    }
  }
  Intens.updateParametersFrom(ParameterValues);
}

} // namespace ComPWA
