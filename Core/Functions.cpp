// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <numeric>
#include "Core/Functions.hpp"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace ComPWA {

bool Inverse::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (checkType != out->type())
    throw BadParameter("Inverse::execute() | Parameter type mismatch!");
  if ((paras.numValues() + paras.numParameters()) != 1)
    throw BadParameter(
        "Inverse::execute() | Expecting a single parameter in list");

  switch (checkType) {
  case ParType::DOUBLE: {
    if (paras.doubleValues().size() == 1) {
      double var = paras.doubleValue(0)->value();
      double result;
      if (var == 0) {
        result = 0;
        LOG(error) << "Inverse::execute() | Division by zero";
      } else
        result = 1 / var;
      out = std::shared_ptr<Parameter>(new Value<double>(out->name(), result));
      break;
    }
    if (paras.doubleParameters().size() == 1) {
      double var = paras.doubleParameter(0)->value();
      double result;
      if (var == 0) {
        result = 0;
        LOG(error) << "Inverse::execute() | Division by zero";
      } else
        result = 1 / var;
      out =
          std::shared_ptr<Parameter>(new DoubleParameter(out->name(), result));
      break;
    }
  }
  default: {
    throw BadParameter("Inverse::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handelt");
  }
  } // end switch
  return true;
} // end execute

bool SquareRoot::execute(ParameterList &paras,
                         std::shared_ptr<Parameter> &out) {
  if (checkType != out->type())
    throw BadParameter("Inverse::SquareRoot() | Parameter type mismatch!");

  if ((paras.numValues() + paras.numParameters()) != 1)
    throw BadParameter(
        "SquareRoot::execute() | Expecting a single parameter in list");

  switch (checkType) {
  case ParType::DOUBLE: {
    if (paras.doubleValues().size() == 1) {
      double var = paras.doubleValue(0)->value();
      out = std::shared_ptr<Parameter>(
          new Value<double>(out->name(), std::sqrt(var)));
      break;
    }
    if (paras.doubleParameters().size() == 1) {
      double var = paras.doubleParameter(0)->value();
      out = std::shared_ptr<Parameter>(
          new DoubleParameter(out->name(), std::sqrt(var)));
      break;
    }
  }
  default: {
    throw BadParameter("SquareRoot::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handelt");
  }
  } // end switch

  return true;
}

struct KahanSummation {
  double sum;
  double correction;
};

/// KahanSummation keeps track of lost bits and reduced the uncertainty in the
/// summation of many large/small numbers.
/// See https://en.wikipedia.org/wiki/Kahan_summation_algorithm
KahanSummation KahanSum(KahanSummation accumulation, double value) {
  KahanSummation result;
  double y = value - accumulation.correction;
  double t = accumulation.sum + y;
  result.correction = (t - accumulation.sum) - y;
  result.sum = t;
  return result;
}

bool AddAll::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (checkType != out->type())
    throw BadParameter("AddAll::SquareRoot() | Parameter type mismatch!");

  switch (checkType) {
  // For multi values we perform a vector addition. All multi values are added
  // and casted to the respective return
  case ParType::MCOMPLEX: {
    size_t n;
    if (paras.mComplexValues().size())
      n = paras.mComplexValue(0)->values().size();
    else if (paras.mDoubleValues().size())
      n = paras.mDoubleValue(0)->values().size();
    else if (paras.mDoubleValues().size())
      n = paras.mIntValue(0)->values().size();
    else
      throw BadParameter(
          "AddAll::execute() | Expecting at least one multi value.");

    // fill MultiComplex parameter
    std::vector<std::complex<double>> results(n, 0.);

    for (auto dv : paras.mIntValues()) {
      if (dv->values().size() != n)
        throw BadParameter("AddAll::execute() | MCOMPLEX: Size of multi int "
                           "value does not match!");
      std::transform(results.begin(), results.end(), dv->values().begin(),
                     results.begin(), std::plus<std::complex<double>>());
    }

    for (auto dv : paras.mDoubleValues()) {
      if (dv->values().size() != n)
        throw BadParameter("AddAll::execute() | MCOMPLEX: Size of multi double "
                           "value does not match!");
      std::transform(results.begin(), results.end(), dv->values().begin(),
                     results.begin(), std::plus<std::complex<double>>());
    }

    for (auto dv : paras.mComplexValues()) {
      if (dv->values().size() != n)
        throw BadParameter(
            "AddAll::execute() | MCOMPLEX: Size of multi complex "
            "value does not match!");
      std::transform(results.begin(), results.end(), dv->values().begin(),
                     results.begin(), std::plus<std::complex<double>>());
    }

    out = std::shared_ptr<Parameter>(
        new Value<std::vector<std::complex<double>>>(out->name(), results));

    break;
  } // end multi complex
  case ParType::MINTEGER: {
    size_t n;
    if (paras.mComplexValues().size())
      throw BadParameter("AddAll::execute() | Return type is double but "
                         "complex value was found!");
    else if (paras.mDoubleValues().size())
      throw BadParameter("AddAll::execute() | Return type is int but "
                         "double value was found!");
    else if (paras.mDoubleValues().size())
      n = paras.mIntValue(0)->values().size();
    else
      throw BadParameter(
          "AddAll::execute() | Expecting at least one multi value.");

    // fill MultiComplex parameter
    std::vector<int> results(n, 0.);

    for (auto dv : paras.mIntValues()) {
      if (dv->values().size() != results.size())
        throw BadParameter("AddAll::execute() | MDOUBLE: Size of multi double "
                           "value does not match!");
      std::transform(results.begin(), results.end(), dv->values().begin(),
                     results.begin(), std::plus<int>());
    }
    out = std::shared_ptr<Parameter>(
        new Value<std::vector<int>>(out->name(), results));
    break;
  } // end multi double
  case ParType::MDOUBLE: {
    size_t n;
    if (paras.mComplexValues().size())
      throw BadParameter("AddAll::execute() | Return type is double but "
                         "complex value was found!");
    else if (paras.mDoubleValues().size())
      n = paras.mDoubleValue(0)->values().size();
    else if (paras.mDoubleValues().size())
      n = paras.mIntValue(0)->values().size();
    else
      throw BadParameter(
          "AddAll::execute() | Expecting at least one multi value.");

    // fill MultiComplex parameter
    std::vector<double> results(n, 0.);

    for (auto dv : paras.mDoubleValues()) {
      if (dv->values().size() != results.size())
        throw BadParameter("AddAll::execute() | MDOUBLE: Size of multi double "
                           "value does not match!");
      std::transform(results.begin(), results.end(), dv->values().begin(),
                     results.begin(), std::plus<double>());
    }
    for (auto dv : paras.mIntValues()) {
      if (dv->values().size() != results.size())
        throw BadParameter("AddAll::execute() | MDOUBLE: Size of multi double "
                           "value does not match!");
      std::transform(results.begin(), results.end(), dv->values().begin(),
                     results.begin(), std::plus<double>());
    }
    out = std::shared_ptr<Parameter>(
        new Value<std::vector<double>>(out->name(), results));
    break;
  } // end multi double
  case ParType::COMPLEX: {
    if (!(paras.complexValues().size() || paras.doubleValues().size() ||
          paras.intValues().size() || paras.boolValues().size()))
      throw BadParameter(
          "AddAll::execute() | COMPLEX: expecting at least one single value!");

    std::complex<double> result(0, 0);
    result = std::accumulate(paras.complexValues().begin(),
                             paras.complexValues().end(), result);
    result = std::accumulate(paras.doubleValues().begin(),
                             paras.doubleValues().end(), result);
    result = std::accumulate(paras.doubleParameters().begin(),
                             paras.doubleParameters().end(), result);
    result = std::accumulate(paras.intValues().begin(), paras.intValues().end(),
                             result);

    // collapse multi values
    for (auto dv : paras.mComplexValues())
      result =
          std::accumulate(dv->values().begin(), dv->values().end(), result);

    for (auto dv : paras.mDoubleValues())
      result =
          std::accumulate(dv->values().begin(), dv->values().end(), result);

    for (auto dv : paras.mIntValues())
      result =
          std::accumulate(dv->values().begin(), dv->values().end(), result);

    out = std::shared_ptr<Parameter>(
        new Value<std::complex<double>>(out->name(), result));
    break;
  } // end complex

  case ParType::DOUBLE: {
    double result = 0;
    result = std::accumulate(paras.doubleValues().begin(),
                             paras.doubleValues().end(), result);
    result = std::accumulate(paras.doubleParameters().begin(),
                             paras.doubleParameters().end(), result);
    result = std::accumulate(paras.intValues().begin(), paras.intValues().end(),
                             result);

    // collapse multi values
    for (auto dv : paras.mDoubleValues()) {
      KahanSummation kaSum = {result};
      auto kaResult = std::accumulate(dv->values().begin(), dv->values().end(),
                                      kaSum, KahanSum);
      result = kaResult.sum;
    }
    for (auto dv : paras.mIntValues()) {
      KahanSummation kaSum = {result};
      auto kaResult = std::accumulate(dv->values().begin(), dv->values().end(),
                                      kaSum, KahanSum);
      result = kaResult.sum;
    }
    out = std::shared_ptr<Parameter>(new Value<double>(out->name(), result));
    break;
  } // end double

  case ParType::INTEGER: {
    double result = 0;
    result = std::accumulate(paras.intValues().begin(), paras.intValues().end(),
                             result);

    // collapse multi values
    for (auto dv : paras.mDoubleValues()) {
      KahanSummation kaSum = {result};
      auto kaResult = std::accumulate(dv->values().begin(), dv->values().end(),
                                      kaSum, KahanSum);
      result = kaResult.sum;
    }
    for (auto dv : paras.mIntValues()) {
      KahanSummation kaSum = {result};
      auto kaResult = std::accumulate(dv->values().begin(), dv->values().end(),
                                      kaSum, KahanSum);
      result = kaResult.sum;
    }

    out = std::shared_ptr<Parameter>(new Value<double>(out->name(), result));
    break;
  } // end int
  default: {
    throw BadParameter("AddAll::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handelt");
  }
  } // end switch
  return true;
}

bool MultAll::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (checkType != out->type())
    throw BadParameter("MultAll::SquareRoot() | Parameter type mismatch!");

  unsigned int nMC = paras.mComplexValues().size();
  unsigned int nMD = paras.mDoubleValues().size();
  unsigned int nMI = paras.mIntValues().size();
  unsigned int nC = paras.complexValues().size();
  unsigned int nD =
      paras.doubleValues().size() + paras.doubleParameters().size();
  unsigned int nI = paras.intValues().size();

  switch (checkType) {

  case ParType::MCOMPLEX: {
    // output multi complex: treat everything non-complex as real,
    // there must be multi complex input
    if (!nMC)
      throw BadParameter("MultAll::execute() | MCOMPLEX: expecting at least "
                         "one multi complex value!");

    std::complex<double> result(1., 0.); // mult up all 1-dim input
    result = std::accumulate(paras.complexValues().begin(),
                             paras.complexValues().end(), result,
                             std::multiplies<std::complex<double>>());
    result = std::accumulate(paras.doubleValues().begin(),
                             paras.doubleValues().end(), result,
                             std::multiplies<std::complex<double>>());
    result = std::accumulate(paras.doubleParameters().begin(),
                             paras.doubleParameters().end(), result,
                             std::multiplies<std::complex<double>>());
    result = std::accumulate(paras.intValues().begin(), paras.intValues().end(),
                             result, std::multiplies<std::complex<double>>());

    // fill MultiComplex parameter
    size_t n = paras.mComplexValue(0)->operator()().size();
    std::vector<std::complex<double>> results(n, result);
    for (auto p : paras.mIntValues()) {
      std::transform(p->operator()().begin(), p->operator()().end(),
                     results.begin(), std::multiplies<std::complex<double>>());
    }
    for (auto p : paras.mDoubleValues()) {
      std::transform(p->operator()().begin(), p->operator()().end(),
                     results.begin(), std::multiplies<std::complex<double>>());
    }
    for (auto p : paras.mComplexValues()) {
      std::transform(p->operator()().begin(), p->operator()().end(),
                     results.begin(), std::multiplies<std::complex<double>>());
    }

    out = std::shared_ptr<Parameter>(
        new Value<std::vector<std::complex<double>>>(out->name(), results));

    break;
  } // end multi complex

  case ParType::MDOUBLE: {
    // output multi double: ignore complex pars, there must be
    // multi double input
    if (!nMD || nMC)
      throw BadParameter(
          "MultAll::execute() | MDOUBLE: Number and/or types do not match");

    double result = 1.;
    result = std::accumulate(paras.doubleValues().begin(),
                             paras.doubleValues().end(), result,
                             std::multiplies<double>());
    result = std::accumulate(paras.doubleParameters().begin(),
                             paras.doubleParameters().end(), result,
                             std::multiplies<double>());
    result = std::accumulate(paras.intValues().begin(), paras.intValues().end(),
                             result, std::multiplies<double>());

    // fill MultiComplex parameter
    size_t n = paras.mDoubleValue(0)->operator()().size();
    std::vector<double> results(n, result);
    for (auto p : paras.mIntValues()) {
      std::transform(p->operator()().begin(), p->operator()().end(),
                     results.begin(), std::multiplies<std::complex<double>>());
    }
    for (auto p : paras.mDoubleValues()) {
      std::transform(p->operator()().begin(), p->operator()().end(),
                     results.begin(), std::multiplies<std::complex<double>>());
    }
    out = std::shared_ptr<Parameter>(
        new Value<std::vector<double>>(out->name(), results));
    break;
  } // end multi double

  case ParType::MINTEGER: {
    // output multi double: ignore complex pars, there must be
    // multi double input
    if (!nMI || nMC || nMD)
      throw BadParameter(
          "MultAll::execute() | MDOUBLE: Number and/or types do not match");

    int result = 1.;
    result = std::accumulate(paras.intValues().begin(), paras.intValues().end(),
                             result, std::multiplies<int>());

    // fill MultiComplex parameter
    size_t n = paras.mIntValue(0)->operator()().size();
    std::vector<int> results(n, result);
    for (auto p : paras.mIntValues()) {
      std::transform(p->operator()().begin(), p->operator()().end(),
                     results.begin(), std::multiplies<int>());
    }
    out = std::shared_ptr<Parameter>(
        new Value<std::vector<int>>(out->name(), results));
    break;
  } // end multi int

  case ParType::COMPLEX: {
    // output complex: collapse everything non-complex as real-part

    if (!nC || nMD || nMC || nMI)
      throw BadParameter("MultAll::execute() | COMPLEX: expecting at least "
                         "one multi complex value!");

    std::complex<double> result(1., 0.); // mult up all 1-dim input
    result = std::accumulate(paras.complexValues().begin(),
                             paras.complexValues().end(), result,
                             std::multiplies<std::complex<double>>());
    result = std::accumulate(paras.doubleValues().begin(),
                             paras.doubleValues().end(), result,
                             std::multiplies<std::complex<double>>());
    result = std::accumulate(paras.doubleParameters().begin(),
                             paras.doubleParameters().end(), result,
                             std::multiplies<std::complex<double>>());
    result = std::accumulate(paras.intValues().begin(), paras.intValues().end(),
                             result, std::multiplies<std::complex<double>>());
    out = std::shared_ptr<Parameter>(
        new Value<std::complex<double>>(out->name(), result));
    break;
  } // end complex

  case ParType::DOUBLE: {
    if (!nD || nC || nMD || nMC || nMI)
      throw BadParameter("MultAll::execute() | DOUBLE: expecting at least "
                         "one multi complex value!");

    double result = 1.; // mult up all 1-dim input
    result = std::accumulate(paras.doubleValues().begin(),
                             paras.doubleValues().end(), result,
                             std::multiplies<double>());
    result = std::accumulate(paras.doubleParameters().begin(),
                             paras.doubleParameters().end(), result,
                             std::multiplies<double>());
    result = std::accumulate(paras.intValues().begin(), paras.intValues().end(),
                             result, std::multiplies<double>());
    out = std::shared_ptr<Parameter>(new Value<double>(out->name(), result));
    break;
  } // end double
  case ParType::INTEGER: {
    if (!nI || nD || nC || nMD || nMC || nMI)
      throw BadParameter("MultAll::execute() | INTEGER: expecting at least "
                         "one multi complex value!");

    int result = 1.; // mult up all 1-dim input
    result = std::accumulate(paras.intValues().begin(), paras.intValues().end(),
                             result, std::multiplies<int>());
    out = std::shared_ptr<Parameter>(new Value<int>(out->name(), result));
    break;
  } // end double
  default: {
    throw BadParameter("MultAll::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handelt");
  }
  } // end switch
  return true;
}

bool LogOf::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (checkType != out->type())
    throw BadParameter("LogOf::execute() | Parameter type mismatch!");

  if (paras.numParameters() != 1)
    throw BadParameter("LogOf::execute() | Expecting only one parameter");

  unsigned int nMC = paras.mComplexValues().size();
  unsigned int nMD = paras.mDoubleValues().size();
  unsigned int nMI = paras.mIntValues().size();
  unsigned int nC = paras.complexValues().size();
  unsigned int nD =
      paras.doubleValues().size() + paras.doubleParameters().size();
  unsigned int nI = paras.intValues().size();

  switch (checkType) {
  case ParType::MDOUBLE: {
    // output multi double: input must be one multi double
    if (!nMD && !nMI)
      throw BadParameter(
          "LogOf::execute() | MDOUBLE: Number and/or types do not match");

    std::vector<double> results;
    std::transform(paras.mDoubleValue(0)->operator()().begin(),
                   paras.mDoubleValue(0)->operator()().end(),
                   std::back_inserter(results),
                   [](double x) { return std::log(x); });
    std::transform(paras.mIntValue(0)->operator()().begin(),
                   paras.mIntValue(0)->operator()().end(),
                   std::back_inserter(results),
                   [](double x) { return std::log(x); });

    out = std::shared_ptr<Parameter>(
        new Value<std::vector<double>>(out->name(), results));
    break;
  } // end multi double

  case ParType::DOUBLE: {
    if (!nD && !nI)
      throw BadParameter(
          "LogOf::execute() | DOUBLE: Number and/or types do not match");

    // output double: log of one double input
    if (paras.doubleValues().size())
      out = std::shared_ptr<Parameter>(new DoubleParameter(
          out->name(), std::log(paras.doubleValue(0)->value())));
    else if (paras.doubleParameters().size())
      out = std::shared_ptr<Parameter>(new DoubleParameter(
          out->name(), std::log(paras.doubleParameter(0)->value())));
    else if (paras.intValues().size())
      out = std::shared_ptr<Parameter>(new DoubleParameter(
          out->name(), std::log(paras.intValue(0)->value())));
    else
      throw std::runtime_error("LogOf::execute() | DOUBLE: something is wrong. "
                               "We should not arrive here!");
    break;
  } // end double
  default: {
    throw BadParameter("LogOf::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handelt");
  }
  } // end switch
  return true;
};

bool Complexify::execute(ParameterList &paras,
                         std::shared_ptr<Parameter> &out) {
  if (checkType != out->type())
    throw BadParameter("Complexify::SquareRoot() | Parameter type mismatch!");

  unsigned int nMC = paras.mComplexValues().size();
  unsigned int nMD = paras.mDoubleValues().size();
  unsigned int nMI = paras.mIntValues().size();
  unsigned int nC = paras.complexValues().size();
  unsigned int nD =
      paras.doubleValues().size() + paras.doubleParameters().size();
  unsigned int nI = paras.intValues().size();

  switch (checkType) {
  case ParType::MCOMPLEX: {
    // output multi complex: input must be two multi double
    if (nMD != 2 || nMC || nMI || nC || nD || nI)
      throw BadParameter(
          "Complexify::execute() | MCOMPLEX: Number and/or types do not match");

    std::vector<std::complex<double>> results;

    // We have to assume here that the magnitude is the first parameter and
    // the phase the second one. We cannot check that.
    std::transform(
        paras.mDoubleValue(0)->operator()().begin(),
        paras.mDoubleValue(0)->operator()().end(),
        paras.mDoubleValue(1)->operator()().begin(), results.begin(),
        [](double r, double phi) { return std::polar(std::abs(r), phi); });

    out = std::shared_ptr<Parameter>(
        new Value<std::vector<std::complex<double>>>(out->name(), results));
    break;
  } // end multi complex
  case ParType::COMPLEX: {
    // output complex: input must be two double
    // output multi complex: input must be two multi double
    if (nD != 2 || nMC || nMD || nMI || nC || nI)
      throw BadParameter(
          "Complexify::execute() | COMPLEX: Number and/or types do not match");

    if (paras.doubleValues().size() == 2) {
        out = std::shared_ptr<Parameter>(new Value<std::complex<double>>(
        out->name(),
        std::polar(std::abs(paras.doubleValue(0)), paras.doubleValue(1)));
    } else if (paras.doubleParameters().size() == 2) {
            out = std::shared_ptr<Parameter>(new Value<std::complex<double>>(
        out->name(),
        std::polar(std::abs(paras.doubleParameter(0)), paras.doubleParameter(1)));
    } else {
      throw std::runtime_error("LogOf::execute() | DOUBLE: something is wrong. "
                               "We should not arrive here!");
    }
    break;
  } // end double
  default: {
    throw BadParameter("Complexify::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handelt");
  }
  } // end switch
  return true;
};

bool ComplexConjugate::execute(ParameterList &paras,
                               std::shared_ptr<Parameter> &out) {
  if (checkType != out->type())
    throw BadParameter(
        "ComplexConjugate::SquareRoot() | Parameter type mismatch!");

  unsigned int nMC = paras.mComplexValues().size();
  unsigned int nMD = paras.mDoubleValues().size();
  unsigned int nMI = paras.mIntValues().size();
  unsigned int nC = paras.complexValues().size();
  unsigned int nD =
      paras.doubleValues().size() + paras.doubleParameters().size();
  unsigned int nI = paras.intValues().size();

  if (nMD || nMI || nD || nI)
    throw BadParameter("ComplexConjugate::execute() | Real numbers given. This "
                       "is mathematically correct but not necessary and "
                       "therefore not implemented.");

  switch (checkType) {
  case ParType::MCOMPLEX: {
    // output complex: input must be one multicomplex
    if (nMC != 1 || nC)
      throw BadParameter("ComplexConjugate::execute() | MCOMPLEX: Number "
                         "and/or types do not match");

    std::vector<std::complex<double>> results;
    std::transform(paras.mComplexValue(0)->operator()().begin(),
                   paras.mComplexValue(0)->operator()().begin(),
                   std::back_inserter(results),
                   [](std::complex<double> c) { return std::conj(c); });

    out = std::shared_ptr<Parameter>(
        new Value<std::vector<std::complex<double>>>(out->name(), results));
    break;
  } // end multi complex
  case ParType::COMPLEX: {
    // output complex: input must be a complex
    if (nC != 1 || nMC)
      throw BadParameter("ComplexConjugate::execute() | COMPLEX: Number and/or "
                         "types do not match");
    out = std::shared_ptr<Parameter>(new Value<std::complex<double>>(
        out->name(), std::conj(paras.complexValue(0)->value())));
    break;
  } // end double
  default: {
    throw BadParameter("ComplexConjugate::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handelt");
  }
  } // end switch

  return true;
};

bool AbsSquare::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (checkType != out->type())
    throw BadParameter("AbsSquare::SquareRoot() | Parameter type mismatch!");

  unsigned int nMC = paras.mComplexValues().size();
  unsigned int nMD = paras.mDoubleValues().size();
  unsigned int nMI = paras.mIntValues().size();
  unsigned int nC = paras.complexValues().size();
  unsigned int nD =
      paras.doubleValues().size() + paras.doubleParameters().size();
  unsigned int nI = paras.intValues().size();

  if (paras.numParameters() != 1)
    throw std::runtime_error("AbsSquare::execute() | Input parameter list "
                             "contains more than one parameter!");

  switch (checkType) {
  case ParType::MINTEGER: {
    if (nMI != 1)
      throw BadParameter("AbsSquare::execute() | MINTEGER: Number and/or "
                         "types do not match");

    size_t n = paras.mIntValue(0)->operator()().size();
    out = std::shared_ptr<Parameter>(
        new Value<std::vector<int>>(out->name(), std::vector<int>(n)));
    auto result = out->operator()();

    std::transform(paras.mDoubleValue(0)->operator()().begin(),
                   paras.mDoubleValue(0)->operator()().begin(), results.begin(),
                   [](int c) { return std::norm(c); });
  }
  case ParType::MDOUBLE: {
    if (nMD == 1) {
      size_t n = paras.mDoubleValue(0)->operator()().size();
      out = std::shared_ptr<Parameter>(
          new Value<std::vector<double>>(out->name(), std::vector<double>(n)));
      auto result = out->operator()();

      std::transform(paras.mDoubleValue(0)->operator()().begin(),
                     paras.mDoubleValue(0)->operator()().begin(),
                     results.begin(), [](double c) { return std::norm(c); });
    } else if (nMC == 1) {
      size_t n = paras.mComplexValue(0)->operator()().size();
      out = std::shared_ptr<Parameter>(
          new Value<std::vector<std::complex<double>>>(
              out->name(), std::vector<std::complex<double>>(n)));
      auto result = out->operator()();

      std::transform(paras.mDoubleValue(0)->operator()().begin(),
                     paras.mDoubleValue(0)->operator()().begin(),
                     results.begin(),
                     [](std::complex<double> c) { return std::norm(c); });
    } else if (nMI == 1) {
      size_t n = paras.mIntValue(0)->operator()().size();
      out = std::shared_ptr<Parameter>(
          new Value<std::vector<int>>(out->name(), std::vector<int>(n)));
      auto result = out->operator()();

      std::transform(paras.mDoubleValue(0)->operator()().begin(),
                     paras.mDoubleValue(0)->operator()().begin(),
                     results.begin(), [](int c) { return std::norm(c); });
    } else {
      throw BadParameter("AbsSquare::execute() | MDOUBLE: Number and/or "
                         "types do not match");
    }
    break;
  } // end multi double
  case ParType::INTEGER: {
    if (nI != 1)
      throw BadParameter("AbsSquare::execute() | INTEGER: Number and/or "
                         "types do not match");
    out = std::shared_ptr<Parameter>(
        new Value<int>(out->name(), std::norm(paras.intValue(0)->value())));
  }
  case ParType::DOUBLE: {
    if (paras.doubleValues().size()) {
      out = std::shared_ptr<Parameter>(new Value<double>(
          out->name(), std::norm(para.doubleValue()->value())));
    } else if (paras.doubleParameters().size()) {
      out = std::shared_ptr<Parameter>(new Value<double>(
          out->name(), std::norm(para.doubleParameters()->value())));
    } else if (nC) {
      out = std::shared_ptr<Parameter>(new Value<double>(
          out->name(), std::norm(para.complexValue()->value())));
    } else {
      throw BadParameter("AbsSquare::execute() | DOUBLE: Number and/or "
                         "types do not match");
    }
    break;
  } // end double
  default: {
    throw BadParameter("AbsSquare::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handelt");
  }
  } // end switch
  return true;
};

bool Power::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (checkType != out->type())
    throw BadParameter("Power::SquareRoot() | Parameter type mismatch!");

  unsigned int nMC = paras.mComplexValues().size();
  unsigned int nMD = paras.mDoubleValues().size();
  unsigned int nMI = paras.mIntValues().size();
  unsigned int nC = paras.complexValues().size();
  unsigned int nD =
      paras.doubleValues().size() + paras.doubleParameters().size();
  unsigned int nI = paras.intValues().size();

  if (nMC + nMD + nD + nI == 0) {
    // TODO: exception no input
    return false;
  }
  // only two double or complex parameter possible
  if (!(nD + nMD + nC + nMC == 2)) {
    // TODO: exception wrong input
    return false;
  }

  switch (checkType) {

  case ParType::MCOMPLEX: {
    // output multi complex: input must be two multi complex
    if (!(nMC == 2)) {
      // TODO: exception wrong input
      return false;
    }

    const std::vector<std::complex<double>> v_base =
        paras.GetMultiComplex(0)->values();
    const std::vector<std::complex<double>> v_exp =
        paras.GetMultiComplex(1)->values();

    std::vector<std::complex<double>> results(v_base.size(),
                                              std::complex<double>(0., 0.));

    for (unsigned int ele = 0; ele < v_base.size(); ele++)
      results.at(ele) = std::pow(v_base.at(ele), v_exp.at(ele));

    out = std::shared_ptr<Parameter>(new MultiComplex(out->name(), results));

    break;
  } // end multi complex

  case ParType::MDOUBLE: {
    // output multi double: input must be two multi double
    if (!(nMD == 2)) {
      // TODO: exception wrong input
      return false;
    }
    // TODO: integer exponent can be calculated much faster
    const std::vector<double> v_base = paras.GetMultiDouble(0)->values();
    const std::vector<double> v_exp = paras.GetMultiDouble(1)->values();

    std::vector<double> results(v_base.size(), 0.);

    for (unsigned int ele = 0; ele < v_base.size(); ele++)
      results.at(ele) = std::pow(v_base.at(ele), v_exp.at(ele));

    out = std::shared_ptr<Parameter>(new MultiDouble(out->name(), results));

    break;
  } // end multi double

  case ParType::COMPLEX: {
    // output complex: power of two complex input
    if (!(nC == 2)) {
      // TODO: exception wrong input
      return false;
    }
    std::shared_ptr<ComplexParameter> tmpA = paras.GetComplexParameter(0);
    std::shared_ptr<ComplexParameter> tmpB = paras.GetComplexParameter(1);
    out = std::shared_ptr<Parameter>(new ComplexParameter(
        out->name(), std::pow(tmpA->value(), tmpB->value())));
    break;
  } // end double

  case ParType::DOUBLE: {
    // output double: power of two double input
    if (!(nD == 2)) {
      // TODO: exception wrong input
      return false;
    }
    std::shared_ptr<DoubleParameter> tmpA = paras.GetDoubleParameter(0);
    std::shared_ptr<DoubleParameter> tmpB = paras.GetDoubleParameter(1);
    out = std::shared_ptr<Parameter>(new DoubleParameter(
        out->name(), std::pow(tmpA->value(), tmpB->value())));
    break;
  } // end double
  default: {
    throw BadParameter("Power::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handelt");
  }
  } // end switch
  return true;
};

} // ns::ComPWA
