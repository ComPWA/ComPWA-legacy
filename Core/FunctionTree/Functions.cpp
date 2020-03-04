// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
#include <functional>
#include <numeric>

#include "Functions.hpp"

#include "ThirdParty/parallelstl/include/pstl/algorithm"
#include "ThirdParty/parallelstl/include/pstl/execution"

namespace ComPWA {
namespace FunctionTree {

void Inverse::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("Inverse::execute() | Parameter type mismatch!");

  if ((paras.doubleValues().size() + paras.doubleParameters().size()) != 1)
    throw BadParameter(
        "Inverse::execute() | Expecting a single parameter or value in list");

  switch (checkType) {
  case ParType::DOUBLE: {
    // Create parameter if not there
    if (!out)
      out = std::make_shared<Value<double>>();
    auto par = std::static_pointer_cast<Value<double>>(out);
    auto &result = par->operator()();
    double var;
    if (paras.doubleValues().size() == 1) {
      var = paras.doubleValue(0)->value();
    } else if (paras.doubleParameters().size() == 1) {
      var = paras.doubleParameter(0)->value();
    } else {
      throw BadParameter(
          "Inverse::execute() | Expecting a single parameter or value in list");
    }
    if (var == 0) {
      result = 0;
      LOG(ERROR) << "Inverse::execute() | Division by zero";
    } else
      result = 1 / var;
    break;
  }
  default: {
    throw BadParameter("Inverse::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handled");
  }
  }
}

void SquareRoot::execute(ParameterList &paras,
                         std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("Inverse::SquareRoot() | Parameter type mismatch!");

  if ((paras.doubleValues().size() + paras.doubleParameters().size()) != 1)
    throw BadParameter(
        "SquareRoot::execute() | Expecting a single parameter in list");

  switch (checkType) {
  case ParType::DOUBLE: {
    // Create parameter if not there
    if (!out)
      out = std::make_shared<Value<double>>();
    auto par = std::static_pointer_cast<Value<double>>(out);
    auto &result = par->operator()();
    double var = 0;
    if (paras.doubleValues().size() == 1) {
      var = paras.doubleValue(0)->value();
    } else if (paras.doubleParameters().size() == 1) {
      var = paras.doubleParameter(0)->value();
    } else {
      throw BadParameter("SquareRoot::execute() | Expecting a single parameter "
                         "or value in list");
    }
    result = std::sqrt(var);
    break;
  }
  default: {
    throw BadParameter("SquareRoot::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handled");
  }
  } // end switch
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

void AddAll::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
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

    if (!out)
      out = MComplex("", n);

    // fill MultiComplex parameter
    auto par =
        std::static_pointer_cast<Value<std::vector<std::complex<double>>>>(out);
    auto &results = par->values(); // reference
    if (results.size() != n) {
      results.resize(n);
    }
    // first get all the scalar inputs and add them
    double initial_real(0.0);
    for (auto const x : paras.doubleValues())
      initial_real += x->value();
    std::complex<double> initial_value(initial_real, 0.0);
    for (auto const x : paras.complexValues())
      initial_value += x->value();
    std::fill(pstl::execution::par_unseq, results.begin(), results.end(),
              initial_value); // reset

    for (auto dv : paras.mComplexValues()) {
      if (dv->values().size() != n)
        throw BadParameter(
            "AddAll::execute() | MCOMPLEX: Size of multi complex "
            "value does not match!");
      std::transform(pstl::execution::par_unseq, results.begin(), results.end(),
                     dv->values().begin(), results.begin(),
                     std::plus<std::complex<double>>());
    }
    for (auto dv : paras.mDoubleValues()) {
      if (dv->values().size() != n)
        throw BadParameter("AddAll::execute() | MCOMPLEX: Size of multi double "
                           "value does not match!");
      std::transform(pstl::execution::par_unseq, results.begin(), results.end(),
                     dv->values().begin(), results.begin(),
                     std::plus<std::complex<double>>());
    }
    for (auto dv : paras.mIntValues()) {
      if (dv->values().size() != n)
        throw BadParameter("AddAll::execute() | MCOMPLEX: Size of multi int "
                           "value does not match!");
      std::transform(pstl::execution::par_unseq, results.begin(), results.end(),
                     dv->values().begin(), results.begin(),
                     std::plus<std::complex<double>>());
    }
    break;
  } // end multi complex
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
    // Create parameter if not there
    if (!out)
      out = MDouble("", n);
    auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
    auto &results = par->values(); // reference
    if (results.size() != n) {
      results.resize(n);
    }
    // first get all the scalar inputs and add them
    double initial_value(0.0);
    for (auto const x : paras.doubleValues())
      initial_value += x->value();
    std::fill(pstl::execution::par_unseq, results.begin(), results.end(),
              initial_value); // reset

    for (auto dv : paras.mDoubleValues()) {
      if (dv->values().size() != results.size())
        throw BadParameter("AddAll::execute() | MDOUBLE: Size of multi double "
                           "value does not match!");
      std::transform(pstl::execution::par_unseq, results.begin(), results.end(),
                     dv->values().begin(), results.begin(),
                     std::plus<double>());
    }
    for (auto dv : paras.mIntValues()) {
      if (dv->values().size() != results.size())
        throw BadParameter("AddAll::execute() | MDOUBLE: Size of multi double "
                           "value does not match!");
      std::transform(pstl::execution::par_unseq, results.begin(), results.end(),
                     dv->values().begin(), results.begin(),
                     std::plus<double>());
    }
    break;
  } // end multi double
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
    // Create parameter if not there
    if (!out)
      out = MInteger("", n);

    auto par = std::static_pointer_cast<Value<std::vector<int>>>(out);
    auto &results = par->values(); // reference
    if (results.size() != n) {
      results.resize(n);
    }
    std::fill(pstl::execution::par_unseq, results.begin(), results.end(),
              0); // reset

    // fill multi integer parameter
    for (auto dv : paras.mIntValues()) {
      if (dv->values().size() != results.size())
        throw BadParameter("AddAll::execute() | MDOUBLE: Size of multi double "
                           "value does not match!");
      std::transform(pstl::execution::par_unseq, results.begin(), results.end(),
                     dv->values().begin(), results.begin(), std::plus<int>());
    }
    break;
  } // end multi double

  case ParType::COMPLEX: {
    if (!(paras.complexValues().size() || paras.doubleValues().size() ||
          paras.intValues().size()))
      throw BadParameter("AddAll::execute() | COMPLEX: expecting at least "
                         "one single value!");
    // Create parameter if not there
    if (!out)
      out = std::make_shared<Value<std::complex<double>>>();
    auto par = std::static_pointer_cast<Value<std::complex<double>>>(out);
    auto &result = par->values();        // reference
    result = std::complex<double>(0, 0); // reset

    for (auto dv : paras.complexValues())
      result += dv->value();
    for (auto dv : paras.doubleValues())
      result += dv->value();
    for (auto dv : paras.doubleParameters())
      result += dv->value();
    for (auto dv : paras.intValues())
      result += dv->value();

    // collapse multi values
    for (auto dv : paras.mComplexValues())
      result +=
          std::accumulate(dv->values().begin(), dv->values().end(), result);

    for (auto dv : paras.mDoubleValues())
      result +=
          std::accumulate(dv->values().begin(), dv->values().end(), result);

    for (auto dv : paras.mIntValues())
      result += std::accumulate(dv->values().begin(), dv->values().end(), 0);

    break;
  } // end complex

  case ParType::DOUBLE: {
    // Create parameter if not there
    if (!out)
      out = std::make_shared<Value<double>>();
    auto par = std::static_pointer_cast<Value<double>>(out);
    auto &result = par->values(); // reference
    result = 0.;                  // reset

    for (auto dv : paras.doubleValues())
      result += dv->value();
    for (auto dv : paras.doubleParameters())
      result += dv->value();
    for (auto dv : paras.intValues())
      result += dv->value();

    // collapse multi values
    for (auto dv : paras.mDoubleValues()) {
      KahanSummation kaSum = {result};
      auto kaResult = std::accumulate(dv->values().begin(), dv->values().end(),
                                      kaSum, KahanSum);
      result += kaResult.sum;
    }
    for (auto dv : paras.mIntValues()) {
      KahanSummation kaSum = {result};
      auto kaResult = std::accumulate(dv->values().begin(), dv->values().end(),
                                      kaSum, KahanSum);
      result += kaResult.sum;
    }
    break;
  } // end double

  case ParType::INTEGER: {
    // Create parameter if not there
    if (!out)
      out = std::make_shared<Value<int>>();
    auto par = std::static_pointer_cast<Value<int>>(out);
    auto &result = par->values(); // reference
    result = 0;                   // reset
    for (auto dv : paras.intValues())
      result += dv->value();

    // collapse multi values
    for (auto dv : paras.mIntValues()) {
      KahanSummation kaSum = {(double)result};
      auto kaResult = std::accumulate(dv->values().begin(), dv->values().end(),
                                      kaSum, KahanSum);
      result += kaResult.sum;
    }
    break;
  } // end int
  default: {
    throw BadParameter("AddAll::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handled");
  }
  } // end switch
}

void MultAll::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("MultAll::execute() | Parameter type mismatch!");

  size_t nMC = paras.mComplexValues().size();
  size_t nMD = paras.mDoubleValues().size();
  size_t nMI = paras.mIntValues().size();
  size_t nC = paras.complexValues().size();
  size_t nD = paras.doubleValues().size() + paras.doubleParameters().size();
  size_t nI = paras.intValues().size();

  switch (checkType) {

  case ParType::MCOMPLEX: {
    // output multi complex: treat everything non-complex as real,
    // there must be multi complex input
    if (!(nMC > 0 || (nMD > 0 && nC > 0)))
      throw BadParameter(
          "MultAll::execute() | MCOMPLEX: expecting at least "
          "one multi complex value or a multi double and a complex scalar!");

    std::complex<double> result(1., 0.); // mult up all 1-dim input
    for (auto p : paras.complexValues())
      result *= p->value();
    for (auto p : paras.doubleValues())
      result *= p->value();
    for (auto p : paras.doubleParameters())
      result *= p->value();
    for (auto p : paras.intValues())
      result *= p->value();

    size_t n(0);
    if (nMC > 0)
      n = paras.mComplexValue(0)->values().size();
    if (n == 0 && nMC == 0) {
      // if there is no multi complex as input but a scalar complex which
      // should be multiplied on a multi double
      n = paras.mDoubleValue(0)->values().size();
    }
    if (!out)
      out = MComplex("", n);
    auto par =
        std::static_pointer_cast<Value<std::vector<std::complex<double>>>>(out);
    auto &results = par->values(); // reference
    if (results.size() != n) {
      results.resize(n);
    }
    std::fill(pstl::execution::par_unseq, results.begin(), results.end(),
              result); // reset

    for (auto p : paras.mComplexValues()) {
      std::transform(pstl::execution::par_unseq, p->values().begin(),
                     p->values().end(), results.begin(), results.begin(),
                     std::multiplies<std::complex<double>>());
    }
    for (auto p : paras.mDoubleValues()) {
      std::transform(pstl::execution::par_unseq, p->values().begin(),
                     p->values().end(), results.begin(), results.begin(),
                     std::multiplies<std::complex<double>>());
    }
    for (auto p : paras.mIntValues()) {
      std::transform(pstl::execution::par_unseq, p->values().begin(),
                     p->values().end(), results.begin(), results.begin(),
                     std::multiplies<std::complex<double>>());
    }
    break;
  } // end multi complex
  case ParType::MDOUBLE: {
    // output multi double: ignore complex pars, there must be
    // multi double input
    if (!nMD || nMC)
      throw BadParameter(
          "MultAll::execute() | MDOUBLE: Number and/or types do not match");

    double result = 1.;
    for (auto p : paras.doubleValues())
      result *= p->value();
    for (auto p : paras.doubleParameters())
      result *= p->value();
    for (auto p : paras.intValues())
      result *= p->value();

    size_t n = paras.mDoubleValue(0)->values().size();
    if (!out)
      out = MDouble("", n);
    // fill MultiComplex parameter
    auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
    auto &results = par->values(); // reference
    if (results.size() != n) {
      results.resize(n);
    }
    std::fill(pstl::execution::par_unseq, results.begin(), results.end(),
              result); // reset

    for (auto p : paras.mDoubleValues()) {
      std::transform(pstl::execution::par_unseq, p->values().begin(),
                     p->values().end(), results.begin(), results.begin(),
                     std::multiplies<double>());
    }
    for (auto p : paras.mIntValues()) {
      std::transform(pstl::execution::par_unseq, p->values().begin(),
                     p->values().end(), results.begin(), results.begin(),
                     std::multiplies<double>());
    }
    break;
  } // end multi double

  case ParType::MINTEGER: {
    // output multi double: ignore complex pars, there must be
    // multi double input
    if (!nMI || nMC || nMD)
      throw BadParameter(
          "MultAll::execute() | MDOUBLE: Number and/or types do not match");

    int result = 1.;
    for (auto p : paras.intValues())
      result *= p->value();

    size_t n = paras.mIntValue(0)->values().size();
    if (!out)
      out = MInteger("", n);

    // fill MultiComplex parameter
    auto par = std::static_pointer_cast<Value<std::vector<int>>>(out);
    auto &results = par->values(); // reference
    if (results.size() != n) {
      results.resize(n);
    }
    std::fill(pstl::execution::par_unseq, results.begin(), results.end(),
              result); // reset

    for (auto p : paras.mIntValues()) {
      std::transform(pstl::execution::par_unseq, p->values().begin(),
                     p->values().end(), results.begin(), results.begin(),
                     std::multiplies<int>());
    }
    break;
  } // end multi int

  case ParType::COMPLEX: {
    // output complex: collapse everything non-complex as real-part

    if (!nC || nMD || nMC || nMI)
      throw BadParameter("MultAll::execute() | COMPLEX: expecting at least "
                         "one multi complex value!");
    if (!out)
      out = std::make_shared<Value<std::complex<double>>>();
    auto par = std::static_pointer_cast<Value<std::complex<double>>>(out);
    auto &result = par->values();          // reference
    result = std::complex<double>(1., 0.); // reset

    for (auto p : paras.complexValues())
      result *= p->value();
    for (auto p : paras.doubleValues())
      result *= p->value();
    for (auto p : paras.doubleParameters())
      result *= p->value();
    for (auto p : paras.intValues())
      result *= p->value();
    break;
  } // end complex

  case ParType::DOUBLE: {
    if (!nD || nC || nMD || nMC || nMI)
      throw BadParameter("MultAll::execute() | DOUBLE: expecting at least "
                         "one multi complex value!");
    if (!out)
      out = std::make_shared<Value<double>>();
    auto par = std::static_pointer_cast<Value<double>>(out);
    auto &result = par->values(); // reference
    result = 1.;                  // reset

    for (auto p : paras.doubleValues())
      result *= p->value();
    for (auto p : paras.doubleParameters())
      result *= p->value();
    for (auto p : paras.intValues())
      result *= p->value();
    break;
  } // end double
  case ParType::INTEGER: {
    if (!nI || nD || nC || nMD || nMC || nMI)
      throw BadParameter("MultAll::execute() | INTEGER: expecting at least "
                         "one multi complex value!");
    if (!out)
      out = std::make_shared<Value<int>>();
    auto par = std::static_pointer_cast<Value<int>>(out);
    auto &result = par->values(); // reference
    result = 1;                   // reset

    for (auto p : paras.intValues())
      result *= p->value();
    break;
  } // end double
  default: {
    throw BadParameter("MultAll::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handled");
  }
  } // end switch
}

void LogOf::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("LogOf::execute() | Parameter type mismatch!");

  if (paras.numParameters() + paras.numValues() != 1)
    throw BadParameter("LogOf::execute() | Expecting only one parameter");

  //  size_t nMC = paras.mComplexValues().size();
  size_t nMD = paras.mDoubleValues().size();
  size_t nMI = paras.mIntValues().size();
  //  size_t nC = paras.complexValues().size();
  size_t nD = paras.doubleValues().size() + paras.doubleParameters().size();
  size_t nI = paras.intValues().size();

  switch (checkType) {
  case ParType::MDOUBLE: {
    // output multi double: input must be one multi double
    if (!nMD && !nMI)
      throw BadParameter(
          "LogOf::execute() | MDOUBLE: Number and/or types do not match");

    if (nMD) {
      size_t n = paras.mDoubleValue(0)->values().size();
      if (!out)
        out = MDouble("", n);
      auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
      auto &results = par->values(); // reference
      if (results.size() != n) {
        results.resize(n);
      }
      std::fill(pstl::execution::par_unseq, results.begin(), results.end(),
                0.); // reset
      std::transform(pstl::execution::par_unseq,
                     paras.mDoubleValue(0)->operator()().begin(),
                     paras.mDoubleValue(0)->operator()().end(), results.begin(),
                     [](double x) { return std::log(x); });
    }
    if (nMI) {
      size_t n = paras.mIntValue(0)->values().size();
      if (!out)
        out = MDouble("", n);
      auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
      auto &results = par->values(); // reference
      if (results.size() != n) {
        results.resize(n);
      }
      std::transform(pstl::execution::par_unseq,
                     paras.mIntValue(0)->operator()().begin(),
                     paras.mIntValue(0)->operator()().end(), results.begin(),
                     [](double x) { return std::log(x); });
    }
    break;
  } // end multi double

  case ParType::DOUBLE: {
    if (!nD && !nI)
      throw BadParameter(
          "LogOf::execute() | DOUBLE: Number and/or types do not match");

    if (!out)
      out = std::make_shared<Value<double>>();
    auto par = std::static_pointer_cast<Value<double>>(out);
    auto &result = par->values(); // reference

    // output double: log of one double input
    if (paras.doubleValues().size())
      result = std::log(paras.doubleValue(0)->value());
    else if (paras.doubleParameters().size())
      result = std::log(paras.doubleParameter(0)->value());
    else if (paras.intValues().size())
      result = std::log(paras.intValue(0)->value());
    else
      throw std::runtime_error("LogOf::execute() | DOUBLE: something is wrong. "
                               "We should not arrive here!");
    break;
  } // end double
  default: {
    throw BadParameter("LogOf::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handled");
  }
  } // end switch
};

void Exp::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("Exp::execute() | Parameter type mismatch!");

  if (paras.numParameters() + paras.numValues() != 1)
    throw BadParameter("Exp::execute() | Expecting only one parameter");

  //  size_t nMC = paras.mComplexValues().size();
  size_t nMD = paras.mDoubleValues().size();
  size_t nMI = paras.mIntValues().size();
  //  size_t nC = paras.complexValues().size();
  size_t nD = paras.doubleValues().size() + paras.doubleParameters().size();
  size_t nI = paras.intValues().size();

  switch (checkType) {
  case ParType::MDOUBLE: {
    // output multi double: input must be one multi double
    if (!nMD && !nMI)
      throw BadParameter(
          "Exp::execute() | MDOUBLE: Number and/or types do not match");

    if (nMD) {
      size_t n = paras.mDoubleValue(0)->values().size();
      if (!out)
        out = MDouble("", n);
      auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
      auto &results = par->values(); // reference
      if (results.size() != n) {
        results.resize(n);
      }
      std::fill(pstl::execution::par_unseq, results.begin(), results.end(),
                0.); // reset
      std::transform(pstl::execution::par_unseq,
                     paras.mDoubleValue(0)->operator()().begin(),
                     paras.mDoubleValue(0)->operator()().end(), results.begin(),
                     [](double x) { return std::exp(x); });
    }
    if (nMI) {
      size_t n = paras.mIntValue(0)->values().size();
      if (!out)
        out = MDouble("", n);
      auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
      auto &results = par->values(); // reference
      if (results.size() != n) {
        results.resize(n);
      }
      std::transform(pstl::execution::par_unseq,
                     paras.mIntValue(0)->operator()().begin(),
                     paras.mIntValue(0)->operator()().end(), results.begin(),
                     [](double x) { return std::exp(x); });
    }
    break;
  } // end multi double

  case ParType::DOUBLE: {
    if (!nD && !nI)
      throw BadParameter(
          "Exp::execute() | DOUBLE: Number and/or types do not match");

    if (!out)
      out = std::make_shared<Value<double>>();
    auto par = std::static_pointer_cast<Value<double>>(out);
    auto &result = par->values(); // reference

    // output double: log of one double input
    if (paras.doubleValues().size())
      result = std::exp(paras.doubleValue(0)->value());
    else if (paras.doubleParameters().size())
      result = std::exp(paras.doubleParameter(0)->value());
    else if (paras.intValues().size())
      result = std::log(paras.intValue(0)->value());
    else
      throw std::runtime_error("Exp::execute() | DOUBLE: something is wrong. "
                               "We should not arrive here!");
    break;
  } // end double
  default: {
    throw BadParameter("Exp::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handled");
  }
  } // end switch
};

void Pow::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("Pow::execute() | Parameter type mismatch!");

  if (paras.numParameters() + paras.numValues() != 1)
    throw BadParameter("Pow::execute() | Expecting only one parameter");

  //  size_t nMC = paras.mComplexValues().size();
  size_t nMD = paras.mDoubleValues().size();
  size_t nMI = paras.mIntValues().size();
  //  size_t nC = paras.complexValues().size();
  size_t nD = paras.doubleValues().size() + paras.doubleParameters().size();
  size_t nI = paras.intValues().size();

  int powerCopy(power);
  switch (checkType) {
  case ParType::MDOUBLE: {
    // output multi double: input must be one multi double
    if (!nMD && !nMI)
      throw BadParameter(
          "Pow::execute() | MDOUBLE: Number and/or types do not match");

    if (nMD) {
      size_t n = paras.mDoubleValue(0)->values().size();
      if (!out)
        out = MDouble("", n);
      auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
      auto &results = par->values(); // reference
      if (results.size() != n) {
        results.resize(n);
      }
      std::fill(pstl::execution::par_unseq, results.begin(), results.end(),
                0.); // reset
      std::transform(pstl::execution::par_unseq,
                     paras.mDoubleValue(0)->operator()().begin(),
                     paras.mDoubleValue(0)->operator()().end(), results.begin(),
                     [powerCopy](double x) { return std::pow(x, powerCopy); });
    }
    if (nMI) {
      size_t n = paras.mIntValue(0)->values().size();
      if (!out)
        out = MDouble("", n);
      auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
      auto &results = par->values(); // reference
      if (results.size() != n) {
        results.resize(n);
      }
      std::transform(pstl::execution::par_unseq,
                     paras.mIntValue(0)->operator()().begin(),
                     paras.mIntValue(0)->operator()().end(), results.begin(),
                     [powerCopy](double x) { return std::pow(x, powerCopy); });
    }
    break;
  } // end multi double

  case ParType::DOUBLE: {
    if (!nD && !nI)
      throw BadParameter(
          "Pow::execute() | DOUBLE: Number and/or types do not match");

    if (!out)
      out = std::make_shared<Value<double>>();
    auto par = std::static_pointer_cast<Value<double>>(out);
    auto &result = par->values(); // reference

    // output double: log of one double input
    if (paras.doubleValues().size())
      result = std::pow(paras.doubleValue(0)->value(), power);
    else if (paras.doubleParameters().size())
      result = std::pow(paras.doubleParameter(0)->value(), power);
    else if (paras.intValues().size())
      result = std::log(paras.intValue(0)->value());
    else
      throw std::runtime_error("Pow::execute() | DOUBLE: something is wrong. "
                               "We should not arrive here!");
    break;
  } // end double
  default: {
    throw BadParameter("Pow::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handled");
  }
  } // end switch
};

void Complexify::execute(ParameterList &paras,
                         std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("Complexify::SquareRoot() | Parameter type mismatch!");

  size_t nMC = paras.mComplexValues().size();
  size_t nMD = paras.mDoubleValues().size();
  size_t nMI = paras.mIntValues().size();
  size_t nC = paras.complexValues().size();
  size_t nD = paras.doubleValues().size() + paras.doubleParameters().size();
  size_t nI = paras.intValues().size();

  switch (checkType) {
  case ParType::MCOMPLEX: {
    // output multi complex: input must be two multi double
    if (nMD != 2 || nMC || nMI || nC || nD || nI)
      throw BadParameter("Complexify::execute() | MCOMPLEX: Number and/or "
                         "types do not match");
    size_t n = paras.mDoubleValue(0)->values().size();
    if (!out)
      out = MComplex("", n);
    auto par =
        std::static_pointer_cast<Value<std::vector<std::complex<double>>>>(out);
    auto &results = par->values(); // reference
    if (results.size() != n) {
      results.resize(n);
    }

    // We have to assume here that the magnitude is the first parameter and
    // the phase the second one. We cannot check that.
    std::transform(
        pstl::execution::par_unseq, paras.mDoubleValue(0)->operator()().begin(),
        paras.mDoubleValue(0)->operator()().end(),
        paras.mDoubleValue(1)->operator()().begin(), results.begin(),
        [](double r, double phi) { return std::polar(std::abs(r), phi); });
    break;
  } // end multi complex
  case ParType::COMPLEX: {
    // output complex: input must be two double
    // output multi complex: input must be two multi double
    if (nD != 2 || nMC || nMD || nMI || nC || nI)
      throw BadParameter("Complexify::execute() | COMPLEX: Number and/or "
                         "types do not match");
    if (!out)
      out = std::make_shared<Value<std::complex<double>>>();
    auto par = std::static_pointer_cast<Value<std::complex<double>>>(out);
    auto &result = par->values(); // reference

    if (paras.doubleValues().size() == 2) {
      result = std::polar(std::abs(paras.doubleValue(0)->value()),
                          paras.doubleValue(1)->value());
    } else if (paras.doubleParameters().size() == 2) {
      result = std::polar(std::abs(paras.doubleParameter(0)->value()),
                          paras.doubleParameter(1)->value());
    } else {
      throw std::runtime_error("LogOf::execute() | DOUBLE: something is wrong. "
                               "We should not arrive here!");
    }
    break;
  } // end double
  default: {
    throw BadParameter("Complexify::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handled");
  }
  } // end switch
};

void ComplexConjugate::execute(ParameterList &paras,
                               std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter(
        "ComplexConjugate::SquareRoot() | Parameter type mismatch!");

  size_t nMC = paras.mComplexValues().size();
  size_t nMD = paras.mDoubleValues().size();
  size_t nMI = paras.mIntValues().size();
  size_t nC = paras.complexValues().size();
  size_t nD = paras.doubleValues().size() + paras.doubleParameters().size();
  size_t nI = paras.intValues().size();

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
    size_t n = paras.mComplexValue(0)->values().size();
    if (!out)
      out = MComplex("", n);
    auto par =
        std::static_pointer_cast<Value<std::vector<std::complex<double>>>>(out);
    auto &results = par->values(); // reference
    if (results.size() != n) {
      results.resize(n);
    }

    std::transform(pstl::execution::par_unseq,
                   paras.mComplexValue(0)->operator()().begin(),
                   paras.mComplexValue(0)->operator()().end(), results.begin(),
                   [](std::complex<double> c) { return std::conj(c); });
    break;
  } // end multi complex
  case ParType::COMPLEX: {
    // output complex: input must be a complex
    if (nC != 1 || nMC)
      throw BadParameter("ComplexConjugate::execute() | COMPLEX: Number and/or "
                         "types do not match");
    if (!out)
      out = std::make_shared<Value<std::complex<double>>>();
    auto par = std::static_pointer_cast<Value<std::complex<double>>>(out);
    auto &result = par->values();                       // reference
    result = std::conj(paras.complexValue(0)->value()); // reset
    break;
  } // end double
  default: {
    throw BadParameter("ComplexConjugate::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handled");
  }
  } // end switch
};

void AbsSquare::execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("AbsSquare::SquareRoot() | Parameter type mismatch!");

  size_t nMC = paras.mComplexValues().size();
  size_t nMD = paras.mDoubleValues().size();
  size_t nMI = paras.mIntValues().size();
  size_t nC = paras.complexValues().size();
  //  size_t nD = paras.doubleValues().size() +
  //  paras.doubleParameters().size();
  size_t nI = paras.intValues().size();

  if (paras.numParameters() + paras.numValues() != 1)
    throw std::runtime_error("AbsSquare::execute() | Input parameter list "
                             "contains more than one parameter!");

  switch (checkType) {

  case ParType::MDOUBLE: {
    if (nMD == 1) {
      size_t n = paras.mDoubleValue(0)->values().size();
      if (!out)
        out = MDouble("", n);
      auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
      auto &results = par->values(); // reference
      if (results.size() != n) {
        results.resize(n);
      }
      std::transform(paras.mDoubleValue(0)->operator()().begin(),
                     paras.mDoubleValue(0)->operator()().end(), results.begin(),
                     [](double c) { return std::norm(c); });
    } else if (nMC == 1) {
      size_t n = paras.mComplexValue(0)->values().size();
      if (!out)
        out = MDouble("", n);
      auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
      auto &results = par->values(); // reference
      if (results.size() != n) {
        results.resize(n);
      }
      std::transform(pstl::execution::par_unseq,
                     paras.mComplexValue(0)->values().begin(),
                     paras.mComplexValue(0)->values().end(), results.begin(),
                     [](std::complex<double> c) { return std::norm(c); });
    } else if (nMI == 1) {
      size_t n = paras.mIntValue(0)->values().size();
      if (!out)
        out = MDouble("", n);
      auto par = std::static_pointer_cast<Value<std::vector<double>>>(out);
      auto &results = par->values(); // reference
      if (results.size() != n) {
        results.resize(n);
      }
      std::transform(pstl::execution::par_unseq,
                     paras.mIntValue(0)->operator()().begin(),
                     paras.mIntValue(0)->operator()().end(), results.begin(),
                     [](int c) { return std::norm(c); });
    } else {
      throw BadParameter("AbsSquare::execute() | MDOUBLE: Number and/or "
                         "types do not match");
    }
    break;
  } // end multi double
  case ParType::MINTEGER: {
    if (nMI != 1)
      throw BadParameter("AbsSquare::execute() | MINTEGER: Number and/or "
                         "types do not match");
    size_t n = paras.mIntValue(0)->values().size();
    if (!out)
      out = MInteger("", n);
    auto par = std::static_pointer_cast<Value<std::vector<int>>>(out);
    auto &results = par->values(); // reference
    if (results.size() != n) {
      results.resize(n);
    }
    std::transform(pstl::execution::par_unseq,
                   paras.mDoubleValue(0)->operator()().begin(),
                   paras.mDoubleValue(0)->operator()().end(), results.begin(),
                   [](int c) { return std::norm(c); });
    break;
  }
  case ParType::INTEGER: {
    if (nI != 1)
      throw BadParameter("AbsSquare::execute() | INTEGER: Number and/or "
                         "types do not match");
    out = std::shared_ptr<Parameter>(
        new Value<int>(out->name(), std::norm(paras.intValue(0)->value())));
    break;
  }
  case ParType::DOUBLE: {
    if (paras.doubleValues().size()) {
      out = std::shared_ptr<Parameter>(new Value<double>(
          out->name(), std::norm(paras.doubleValue(0)->value())));
    } else if (paras.doubleParameters().size()) {
      out = std::shared_ptr<Parameter>(new Value<double>(
          out->name(), std::norm(paras.doubleParameter(0)->value())));
    } else if (nC) {
      out = std::shared_ptr<Parameter>(new Value<double>(
          out->name(), std::norm(paras.complexValue(0)->value())));
    } else {
      throw BadParameter("AbsSquare::execute() | DOUBLE: Number and/or "
                         "types do not match");
    }
    break;
  } // end double
  default: {
    throw BadParameter("AbsSquare::execute() | Parameter of type " +
                       std::to_string(checkType) + " can not be handled");
  }
  } // end switch
};

} // namespace FunctionTree
} // namespace ComPWA
