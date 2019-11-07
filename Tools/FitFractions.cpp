#include "FitFractions.hpp"

#include "Core/FitResult.hpp"
#include "Core/Function.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Data/DataSet.hpp"
#include "Tools/Integration.hpp"

#include <map>

namespace ComPWA {
namespace Tools {

std::ostream &operator<<(std::ostream &os, const FitFractionList &FFList) {
  os << "\nFit Fractions:\n";
  for (auto const &ff : FFList) {
    os << ff.Name << ": " << ff.Value << " +- " << ff.Error << "\n";
  }
  return os;
}

/// Calculates the fit fractions with errors via error propagation from the
/// covariance matrix. The gradients are calculated via numerical
/// differentiation:
/// \f[
/// f´(x) = \frac{f(x+h) - f(x-h)}{2h} + O(h^2)
/// \f]
FitFractionList
FitFractions::calculateFitFractionsWithCovarianceErrorPropagation(
    const std::vector<std::pair<IntensityComponent, IntensityComponent>>
        &Components,
    const ComPWA::Data::DataSet &PhspSample, const ComPWA::FitResult &Result) {

  FitFractionList FitFractions;

  for (auto Component : Components) {
    auto NumeratorData = getIntegralData(Component.first, PhspSample, Result);
    auto DenominatorData =
        getIntegralData(Component.second, PhspSample, Result);

    double FitFractionValue =
        std::get<1>(NumeratorData) / std::get<1>(DenominatorData);

    double FitFractionError = 0.;
    if (Result.IsValid) {
      // FF = N/D -> calculate Jacobi matrix (because we have a single ff, this
      // is just a vector)
      std::vector<double> JacobiColumn;
      std::vector<std::vector<double>> CovMatrix;
      std::tie(JacobiColumn, CovMatrix) = buildJacobiAndCovariance(
          std::make_tuple(std::get<1>(NumeratorData),
                          std::get<2>(NumeratorData)),
          std::make_tuple(std::get<1>(DenominatorData),
                          std::get<2>(DenominatorData)),
          Result);

      // calculate J^T x Cov x J
      for (unsigned int i = 0; i < JacobiColumn.size(); ++i) {
        for (unsigned int j = 0; j < JacobiColumn.size(); ++j) {
          FitFractionError +=
              Result.CovarianceMatrix[i][j] * JacobiColumn[i] * JacobiColumn[j];
        }
      }
    } else {
      LOG(INFO) << "FitFractions::"
                   "calculateFitFractionsWithCovarianceErrorPropagation() | No "
                   "valid fit result. Skip error calculation.";
    }

    FitFractions.push_back(
        {std::get<0>(NumeratorData) + "/" + std::get<0>(DenominatorData),
         FitFractionValue, FitFractionError});
  }

  // the fit fraction errors are currently the squared values, take a sqrt
  for (auto &f : FitFractions) {
    f.Error = std::sqrt(f.Error);
    LOG(TRACE) << "calculateFitFractionsWithCovarianceErrorPropagation(): fit "
                  "fraction for ("
               << f.Name << ") is " << f.Value << " +- " << f.Error;
  }

  return FitFractions;
}

std::tuple<std::vector<double>, std::vector<std::vector<double>>>
FitFractions::buildJacobiAndCovariance(
    const std::tuple<double, std::vector<DerivativeData>> &NominatorDerivatives,
    const std::tuple<double, std::vector<DerivativeData>>
        &DenominatorDerivatives,
    const ComPWA::FitResult &Result) {
  std::vector<double> Column;
  std::vector<size_t> ParameterPositions;

  double N = std::get<0>(NominatorDerivatives);
  double D = std::get<0>(DenominatorDerivatives);

  auto dN = std::get<1>(NominatorDerivatives);
  auto dD = std::get<1>(DenominatorDerivatives);

  size_t PositionCounter(0);
  for (auto x : Result.FinalParameters) {
    if (x.IsFixed)
      continue;
    auto FoundNominator =
        std::find_if(dN.begin(), dN.end(), [&x](auto const &par) {
          return par.ParameterName == x.Name;
        });
    auto FoundDenominator =
        std::find_if(dD.begin(), dD.end(), [&x](auto const &par) {
          return par.ParameterName == x.Name;
        });
    if (dN.end() != FoundNominator && dD.end() != FoundDenominator) {
      // if both derivatives != 0
      Column.push_back((FoundNominator->ValueAtParameterPlusEpsilon /
                            FoundDenominator->ValueAtParameterPlusEpsilon -
                        FoundNominator->ValueAtParameterMinusEpsilon /
                            FoundDenominator->ValueAtParameterMinusEpsilon) /
                       FoundNominator->StepSize);
    } else if (dN.end() != FoundNominator) {
      // if only nominator derivative != 0
      Column.push_back((FoundNominator->ValueAtParameterPlusEpsilon / D -
                        FoundNominator->ValueAtParameterMinusEpsilon / D) /
                       FoundNominator->StepSize);
    } else if (dD.end() != FoundDenominator) {
      // if only denominator derivative != 0
      Column.push_back((N / FoundDenominator->ValueAtParameterPlusEpsilon -
                        N / FoundDenominator->ValueAtParameterMinusEpsilon) /
                       FoundDenominator->StepSize);
    } else {
      ++PositionCounter;
      continue;
    }
    ParameterPositions.push_back(PositionCounter);
    ++PositionCounter;
  }

  std::vector<std::vector<double>> SubCovarianceMatrix(
      ParameterPositions.size(),
      std::vector<double>(ParameterPositions.size()));
  for (size_t i = 0; i < ParameterPositions.size(); ++i) {
    for (size_t j = 0; j < ParameterPositions.size(); ++j) {
      SubCovarianceMatrix[i][j] =
          Result.CovarianceMatrix[ParameterPositions[i]][ParameterPositions[j]];
    }
  }

  return std::make_tuple(Column, SubCovarianceMatrix);
}

std::tuple<std::string, double, std::vector<FitFractions::DerivativeData>>
FitFractions::getIntegralData(IntensityComponent IntensComponent,
                              const ComPWA::Data::DataSet &PhspSample,
                              const ComPWA::FitResult &Result) {
  std::string Name(IntensComponent.first);

  // only calculated the integrals if they have not once before
  auto FoundIntegrals = IntensityGradientDataMapping.find(Name);

  if (IntensityGradientDataMapping.end() == FoundIntegrals) {
    ComPWA::initializeWithFitResult(*IntensComponent.second, Result);

    FoundIntegrals =
        IntensityGradientDataMapping
            .insert({Name, calculateIntensityIntegralData(
                               *IntensComponent.second, PhspSample, Result)})
            .first;
  }
  return std::tuple_cat(std::make_tuple(FoundIntegrals->first),
                        FoundIntegrals->second);
}

std::tuple<double, std::vector<FitFractions::DerivativeData>>
FitFractions::calculateIntensityIntegralData(
    ComPWA::Intensity &Intens, const ComPWA::Data::DataSet &Sample,
    const ComPWA::FitResult &Result) {

  double Integral = ComPWA::Tools::integrate(Intens, Sample);

  std::vector<Parameter> TempParameters = Intens.getParameters();
  std::vector<double> TempParameterValues;
  for (auto const &x : TempParameters)
    TempParameterValues.push_back(x.Value);

  std::vector<DerivativeData> Derivatives;
  for (unsigned int i = 0; i < TempParameters.size(); ++i) {
    auto found = std::find_if(Result.FinalParameters.begin(),
                              Result.FinalParameters.end(),
                              [refname = TempParameters[i].Name](
                                  auto const &p) { return p.Name == refname; });
    if (Result.FinalParameters.end() != found) {
      if (found->IsFixed) {
        continue;
      }
      double TempValue = found->Value;

      // a step size of sqrt(epsilon) * x produces small rounding errors
      double h = std::sqrt(std::numeric_limits<double>::epsilon());
      if (TempValue != 0.0)
        h *= TempValue;

      DerivativeData GradientData;
      GradientData.ParameterName = TempParameters[i].Name;
      double up(TempValue + h);
      TempParameterValues[i] = up;
      Intens.updateParametersFrom(TempParameterValues);
      GradientData.ValueAtParameterPlusEpsilon =
          ComPWA::Tools::integrate(Intens, Sample);

      double down(TempValue - h);
      TempParameterValues[i] = down;
      GradientData.StepSize = up - down;
      // if the compiler optimizes the above line with math associativity
      // (by setting it to 2*h), then we get an additional numerical error
      Intens.updateParametersFrom(TempParameterValues);
      GradientData.ValueAtParameterMinusEpsilon =
          ComPWA::Tools::integrate(Intens, Sample);

      Derivatives.push_back(GradientData);

      // reset values
      TempParameterValues[i] = TempValue;
    }
  }
  return std::make_tuple(Integral, Derivatives);
}

} // namespace Tools
} // namespace ComPWA
