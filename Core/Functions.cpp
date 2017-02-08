/*
 * Functions.cpp
 *
 *  Created on: 8 Sep 2016
 *      Author: weidenka
 */

#include "Core/Functions.hpp"

namespace ComPWA {

bool Inverse::execute(ParameterList &paras,
                      std::shared_ptr<AbsParameter> &out) {
  if (checkType != out->type())
    return false;
  unsigned int nMC = paras.GetNMultiComplex();
  unsigned int nMD = paras.GetNMultiDouble();
  unsigned int nC = paras.GetNComplex();
  unsigned int nD = paras.GetNDouble();
  unsigned int nI = paras.GetNInteger();

  if (nMC + nC != 0) {
    // TODO: can't handle complex numbers
    return false;
  }
  if (nMD + nD + nI != 1) {
    // TODO: exception too many variables
    return false;
  }

  switch (checkType) {
  case ParType::DOUBLE: {
    double var = paras.GetDoubleParameterValue(0);
    if (var == 0) {
      out =
          std::shared_ptr<AbsParameter>(new DoubleParameter(out->GetName(), 0));
      // TODO: exception dividion by 0
      LOG(error) << "Inverse::execute() | Division by zero";
    } else
      out = std::shared_ptr<AbsParameter>(
          new DoubleParameter(out->GetName(), 1 / var));
    break;
  } // end double
  default: {
    // TODO: exception output partype wrong
    return false;
  }
  } // end switch
  return true;
} // end execute

bool SquareRoot::execute(ParameterList &paras,
                         std::shared_ptr<AbsParameter> &out) {
  if (checkType != out->type())
    return false;
  unsigned int nMC = paras.GetNMultiComplex();
  unsigned int nMD = paras.GetNMultiDouble();
  unsigned int nC = paras.GetNComplex();
  unsigned int nD = paras.GetNDouble();
  unsigned int nI = paras.GetNInteger();

  if (nMC + nC != 0) {
    // TODO: can't handle complex numbers
    return false;
  }
  if (nMD + nD + nI != 1) {
    // TODO: exception too many variables
    return false;
  }

  switch (checkType) {
  case ParType::DOUBLE: {
    double var = paras.GetDoubleParameterValue(0);
    if (var < 0) {
      out = std::shared_ptr<AbsParameter>(
          new DoubleParameter(out->GetName(), -1));
      // TODO: exception argument <0
      LOG(error) << "SquareRoot::execute() | Argument "
                    "negative! Returning -1";
    } else
      out = std::shared_ptr<AbsParameter>(
          new DoubleParameter(out->GetName(), sqrt(var)));
    break;
  } // end double
  default: {
    // TODO: exception output partype wrong
    return false;
  }
  } // end switch
  return true;
}

bool AddAll::execute(ParameterList &paras, std::shared_ptr<AbsParameter> &out) {
  if (checkType != out->type())
    return false;
  unsigned int nMC = paras.GetNMultiComplex();
  unsigned int nMD = paras.GetNMultiDouble();
  unsigned int nC = paras.GetNComplex();
  unsigned int nD = paras.GetNDouble();
  unsigned int nI = paras.GetNInteger();

  if (nMC + nMD + nD + nI + nC == 0) {
    // TODO: exception no input
    return false;
  }

  switch (checkType) {

  case ParType::MCOMPLEX: {
    // output multi complex: treat everything non-complex as real,
    // there must be multi complex input
    if (!nMC) {
      // TODO: exception wrong input
      return false;
    }

    unsigned int nElements = paras.GetMultiComplex(0)->GetNValues();

    std::complex<double> result(0, 0); // sum up all 1-dim input
    // sum up complex parameter
    for (unsigned int i = 0; i < nC; i++) {
      result += paras.GetComplexParameter(i)->GetValue();
    }
    // sum up double parameter
    for (unsigned int i = 0; i < nD; i++) {
      result += paras.GetDoubleParameter(i)->GetValue();
    }
    // sum up integer parameter
    for (unsigned int i = 0; i < nI; i++) {
      result += paras.GetIntegerParameter(i)->GetValue();
    }

    // fill MultiComplex parameter
    std::vector<std::complex<double>> results(nElements, result);
    for (unsigned int i = 0; i < nMD; i++) {
      const std::vector<double> v_tmp = paras.GetMultiDouble(i)->GetValues();
      for (unsigned int ele = 0; ele < v_tmp.size(); ele++)
        results[ele] += v_tmp.at(ele);
    }
    for (unsigned int i = 0; i < nMC; i++) {
      const std::vector<std::complex<double>> v_tmp =
          paras.GetMultiComplex(i)->GetValues();
      for (unsigned int ele = 0; ele < v_tmp.size(); ele++)
        results[ele] += v_tmp.at(ele);
    }

    out = std::shared_ptr<AbsParameter>(
        new MultiComplex(out->GetName(), results));

    break;
  } // end multi complex

  case ParType::MDOUBLE: {
    // output multi double: ignore complex pars, there must be
    // multi double input
    if (!nMD) {
      // TODO: exception wrong input
      return false;
    }
    unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
    double result = 0; // sum up all 1-dim input
    // sum up double parameter
    for (unsigned int i = 0; i < nD; i++) {
      result += paras.GetDoubleParameter(i)->GetValue();
    }
    // sum up integer parameter
    for (unsigned int i = 0; i < nI; i++) {
      result += paras.GetIntegerParameter(i)->GetValue();
    }

    // fill MultiDouble parameter
    std::vector<double> results(nElements, result);
    for (unsigned int i = 0; i < nMD; i++) {
      std::vector<double> v_tmp = paras.GetMultiDouble(i)->GetValues();
      for (unsigned int ele = 0; ele < v_tmp.size(); ele++)
        results[ele] += v_tmp.at(ele);
    }
    out =
        std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(), results));

    break;
  } // end multi double

  case ParType::COMPLEX: {
    // output complex: collapse everything non-complex as real-part
    std::complex<double> result(0, 0);

    // sum up complex parameter
    for (unsigned int i = 0; i < nC; i++) {
      result += paras.GetComplexParameter(i)->GetValue();
    }
    // sum up double parameter
    for (unsigned int i = 0; i < nD; i++) {
      result += paras.GetDoubleParameter(i)->GetValue();
    }
    // sum up integer parameter
    for (unsigned int i = 0; i < nI; i++) {
      result += paras.GetIntegerParameter(i)->GetValue();
    }
    // collapse MultiComplex parameter
    for (unsigned int i = 0; i < nMC; i++) {
      std::shared_ptr<MultiComplex> tmp = paras.GetMultiComplex(i);
      for (unsigned int ele = 0; ele < tmp->GetNValues(); ele++)
        result += tmp->GetValue(ele);
    }
    // collapse MultiDoubles parameter
    for (unsigned int i = 0; i < nMD; i++) {
      std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
      for (unsigned int ele = 0; ele < tmp->GetNValues(); ele++)
        result += tmp->GetValue(ele);
    }

    out = std::shared_ptr<AbsParameter>(
        new ComplexParameter(out->GetName(), result));
    break;
  } // end complex

  case ParType::DOUBLE: {
    // output double: ignore complex pars, collapse everything else
    double result = 0;

    // sum up double parameter
    for (unsigned int i = 0; i < nD; i++) {
      result += paras.GetDoubleParameter(i)->GetValue();
    }
    // sum up integer parameter
    for (unsigned int i = 0; i < nI; i++) {
      result += paras.GetIntegerParameter(i)->GetValue();
    }
    // collapse MultiDoubles parameter
    for (unsigned int i = 0; i < nMD; i++) {
      std::shared_ptr<MultiDouble> tmp = paras.GetMultiDouble(i);
      for (unsigned int ele = 0; ele < tmp->GetNValues(); ele++)
        result += tmp->GetValue(ele);
    }

    out = std::shared_ptr<AbsParameter>(
        new DoubleParameter(out->GetName(), result));
    break;
  } // end double

  default: {
    // TODO: exception output partype wrong
    return false;
  }

  } // end switch
  return true;
}

bool MultAll::execute(ParameterList &paras,
                      std::shared_ptr<AbsParameter> &out) {
  if (checkType != out->type())
    return false;

  unsigned int nMC = paras.GetNMultiComplex();
  unsigned int nMD = paras.GetNMultiDouble();
  unsigned int nC = paras.GetNComplex();
  unsigned int nD = paras.GetNDouble();
  unsigned int nI = paras.GetNInteger();

  if (nMC + nMD + nD + nI + nC == 0) {
    // TODO: exception no input
    return false;
  }

  switch (checkType) {

  case ParType::MCOMPLEX: {
    // output multi complex: treat everything non-complex as real,
    // there must be multi complex input
    if (!nMC) {
      // TODO: exception wrong input
      return false;
    }

    unsigned int nElements = paras.GetMultiComplex(0)->GetNValues();

    std::complex<double> result(1., 0.); // mult up all 1-dim input
    // mult up complex parameter
    for (unsigned int i = 0; i < nC; i++) {
      result *= paras.GetComplexParameter(i)->GetValue();
    }
    // mult up double parameter
    for (unsigned int i = 0; i < nD; i++) {
      result *= paras.GetDoubleParameter(i)->GetValue();
    }
    // mult up integer parameter
    for (unsigned int i = 0; i < nI; i++) {
      result *= paras.GetIntegerParameter(i)->GetValue();
    }

    // fill MultiComplex parameter
    std::vector<std::complex<double>> results(nElements, result);
    for (unsigned int i = 0; i < nMD; i++) {
      const std::vector<double> tmpVec = paras.GetMultiDouble(i)->GetValues();
      for (unsigned int ele = 0; ele < tmpVec.size(); ele++)
        results[ele] *= tmpVec.at(ele);
    }
    for (unsigned int i = 0; i < nMC; i++) {
      const std::vector<std::complex<double>> v_tmp =
          paras.GetMultiComplex(i)->GetValues();
      for (unsigned int ele = 0; ele < v_tmp.size(); ele++)
        results[ele] *= v_tmp.at(ele);
    }

    out = std::shared_ptr<AbsParameter>(
        new MultiComplex(out->GetName(), results));

    break;
  } // end multi complex

  case ParType::MDOUBLE: {
    // output multi double: ignore complex pars, there must be
    // multi double input
    if (!nMD) {
      // TODO: exception wrong input
      return false;
    }
    unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
    double result = 1.; // sum up all 1-dim input
    // mult up double parameter
    for (unsigned int i = 0; i < nD; i++) {
      result *= paras.GetDoubleParameter(i)->GetValue();
    }
    // mult up integer parameter
    for (unsigned int i = 0; i < nI; i++) {
      result *= paras.GetIntegerParameter(i)->GetValue();
    }
    // fill MultiDouble parameter
    std::vector<double> results(nElements, result);
    for (unsigned int i = 0; i < nMD; i++) {
      std::vector<double> v_tmp = paras.GetMultiDouble(i)->GetValues();
      for (unsigned int ele = 0; ele < v_tmp.size(); ele++)
        results[ele] *= v_tmp.at(ele);
    }
    out =
        std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(), results));

    break;
  } // end multi double

  case ParType::COMPLEX: {
    // output complex: collapse everything non-complex as real-part
    std::complex<double> result(1., 0);

    // mult up complex parameter
    for (unsigned int i = 0; i < nC; i++) {
      result *= paras.GetComplexParameter(i)->GetValue();
    }
    // mult up double parameter
    for (unsigned int i = 0; i < nD; i++) {
      result *= paras.GetDoubleParameter(i)->GetValue();
    }
    // mult up integer parameter
    for (unsigned int i = 0; i < nI; i++) {
      result *= paras.GetIntegerParameter(i)->GetValue();
    }
    // collapse MultiComplex parameter
    for (unsigned int i = 0; i < nMC; i++) {
      std::vector<std::complex<double>> v_tmp =
          paras.GetMultiComplex(i)->GetValues();
      for (unsigned int ele = 0; ele < v_tmp.size(); ele++)
        result *= v_tmp.at(ele);
    }
    // collapse MultiDoubles parameter
    for (unsigned int i = 0; i < nMD; i++) {
      std::vector<double> v_tmp = paras.GetMultiDouble(i)->GetValues();
      for (unsigned int ele = 0; ele < v_tmp.size(); ele++)
        result *= v_tmp.at(ele);
    }

    out = std::shared_ptr<AbsParameter>(
        new ComplexParameter(out->GetName(), result));
    break;
  } // end complex

  case ParType::DOUBLE: {
    // output double: ignore complex pars, collapse everything else
    double result = 1.;

    // mult up double parameter
    for (unsigned int i = 0; i < nD; i++) {
      result *= paras.GetDoubleParameter(i)->GetValue();
    }
    // mult up integer parameter
    for (unsigned int i = 0; i < nI; i++) {
      result *= paras.GetIntegerParameter(i)->GetValue();
    }
    // collapse MultiDoubles parameter
    for (unsigned int i = 0; i < nMD; i++) {
      std::vector<double> v_tmp = paras.GetMultiDouble(i)->GetValues();
      for (unsigned int ele = 0; ele < v_tmp.size(); ele++)
        result *= v_tmp.at(ele);
    }

    out = std::shared_ptr<AbsParameter>(
        new DoubleParameter(out->GetName(), result));
    break;
  } // end double

  default: {
    // TODO: exception output partype wrong
    return false;
  }

  } // end switch

  return true;
}

bool LogOf::execute(ParameterList &paras, std::shared_ptr<AbsParameter> &out) {

  if (checkType != out->type())
    return false;

  //    unsigned int nMC = paras.GetNMultiComplex();
  unsigned int nMD = paras.GetNMultiDouble();
  //    unsigned int nC = paras.GetNComplex();
  unsigned int nD = paras.GetNDouble();
  //    unsigned int nI = paras.GetNInteger();

  if (nMD + nD == 0) {
    // TODO: exception no input
    return false;
  }
  // only one parameter possible
  if ((nMD + nD) > 1) {
    // TODO: exception wrong input
    return false;
  }

  switch (checkType) {

  case ParType::MDOUBLE: {
    // output multi double: input must be one multi double
    if (!nMD) {
      // TODO: exception wrong input
      return false;
    }
    // fill MultiDouble parameter
    const std::vector<double> tmp = paras.GetMultiDouble(0)->GetValues();
    std::vector<double> results(tmp.size(), 0.);
    for (unsigned int ele = 0; ele < tmp.size(); ele++)
      results[ele] = std::log(tmp.at(ele));

    out =
        std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(), results));
    break;
  } // end multi double

  case ParType::DOUBLE: {
    if (!nD) {
      // TODO: exception wrong input
      return false;
    }
    // output double: log of one double input
    out = std::shared_ptr<AbsParameter>(new DoubleParameter(
        out->GetName(), std::log(paras.GetDoubleParameterValue(0))));
    break;
  } // end double

  default: {
    // TODO: exception output partype wrong
    return false;
  }

  } // end switch

  return true;
};

bool Complexify::execute(ParameterList &paras,
                         std::shared_ptr<AbsParameter> &out) {
  if (checkType != out->type())
    return false;

  unsigned int nMC = paras.GetNMultiComplex();
  unsigned int nMD = paras.GetNMultiDouble();
  unsigned int nC = paras.GetNComplex();
  unsigned int nD = paras.GetNDouble();
  unsigned int nI = paras.GetNInteger();

  if (nMD + nD == 0) {
    // TODO: exception no input
    return false;
  }
  // only one parameter possible
  if ((nMD + nD) != 2 || (nMC + nC + nI) != 0) {
    // TODO: exception wrong input
    return false;
  }

  switch (checkType) {

  case ParType::MCOMPLEX: {
    // output multi complex: input must be two multi double
    if (!(nMD == 2)) {
      // TODO: exception wrong input
      return false;
    }
    // fill MultiDouble parameter
    const std::vector<double> v_mag = paras.GetMultiDouble(0)->GetValues();
    const std::vector<double> v_phi = paras.GetMultiDouble(1)->GetValues();
    if (v_mag.size() == v_phi.size())
      throw std::runtime_error("Complexify::execute() | "
                               "Vector sizes do not match!");

    std::vector<std::complex<double>> results(v_mag.size(),
                                              std::complex<double>(0., 0.));

    for (unsigned int ele = 0; ele < v_mag.size(); ele++)
      results[ele] = std::complex<double>(
          std::fabs(v_mag.at(ele)) * std::cos(v_phi.at(ele)), // a*cos(phi)
          std::fabs(v_mag.at(ele)) * std::sin(v_phi.at(ele))  // a*sin(phi)
          );

    out = std::shared_ptr<AbsParameter>(
        new MultiComplex(out->GetName(), results));

    break;
  } // end multi complex

  case ParType::COMPLEX: {
    // output complex: input must be two double
    if (!(nD == 2)) {
      // TODO: exception wrong input
      return false;
    }
    double a = std::fabs(paras.GetDoubleParameter(0)->GetValue());
    double phi = paras.GetDoubleParameter(1)->GetValue();
    out = std::shared_ptr<AbsParameter>(new ComplexParameter(
        out->GetName(),
        std::complex<double>(a * std::cos(phi), a * std::sin(phi))));
    break;
  } // end double

  default: {
    // TODO: exception output partype wrong
    return false;
  }

  } // end switch

  return true;
};

bool AbsSquare::execute(ParameterList &paras,
                        std::shared_ptr<AbsParameter> &out) {
  if (checkType != out->type())
    return false;

  unsigned int nMC = paras.GetNMultiComplex();
  unsigned int nMD = paras.GetNMultiDouble();
  unsigned int nC = paras.GetNComplex();
  unsigned int nD = paras.GetNDouble();
  unsigned int nI = paras.GetNInteger();

  if (nMC + nMD + nD + nI + nC == 0)
    throw std::runtime_error("AbsSquare::execute() | Input parameter list "
                             "is empty!");

  // only one parameter possible
  if ((nMC + nMD + nD + nI + nC) > 1)
    throw std::runtime_error("AbsSquare::execute() | More then one type "
                             "of parameter in input!");

  switch (checkType) {
  case ParType::MDOUBLE: {
    // output multi double: input must be multi double or multi complex
    if (!nMD && !nMC)
      throw std::runtime_error(
          "AbsSquare::execute() | Requested output is"
          " MDOUBLE but no MDOUBLE or MCOMPLEX was given as input!");
    if (nMD) {
      const std::vector<double> v_tmp = paras.GetMultiDouble(0)->GetValues();

      std::vector<double> results(v_tmp.size(), 0.);
      for (unsigned int ele = 0; ele < v_tmp.size(); ele++)
        results.at(ele) = std::norm(v_tmp.at(ele));

      out = std::shared_ptr<AbsParameter>(
          new MultiDouble(out->GetName(), results));
    } else if (nMC) {
      const std::vector<std::complex<double>> v_tmp =
          paras.GetMultiComplex(0)->GetValues();

      std::vector<double> results(v_tmp.size(), 0.);
      for (unsigned int ele = 0; ele < v_tmp.size(); ele++)
        results.at(ele) = std::norm(v_tmp.at(ele));

      out = std::shared_ptr<AbsParameter>(
          new MultiDouble(out->GetName(), results));
    }

    break;
  } // end multi double

  case ParType::DOUBLE: {
    // output double: norm of one double or complex input
    if (!nD && !nC)
      throw std::runtime_error(
          "AbsSquare::execute() | Requested output "
          "is DOUBLE but no DOUBLE or COMPLEX was given as input!");
    if (nD) {
      std::shared_ptr<DoubleParameter> tmp = paras.GetDoubleParameter(0);
      out = std::shared_ptr<AbsParameter>(
          new DoubleParameter(out->GetName(), std::norm(tmp->GetValue())));
    } else if (nC) {
      std::shared_ptr<ComplexParameter> tmp = paras.GetComplexParameter(0);
      out = std::shared_ptr<AbsParameter>(
          new DoubleParameter(out->GetName(), std::norm(tmp->GetValue())));
    }
    break;
  } // end double

  default: {
    // TODO: exception output partype wrong
    return false;
  }

  } // end switch

  return true;
};

bool Power::execute(ParameterList &paras, std::shared_ptr<AbsParameter> &out) {
  if (checkType != out->type())
    return false;

  unsigned int nMC = paras.GetNMultiComplex();
  unsigned int nMD = paras.GetNMultiDouble();
  unsigned int nC = paras.GetNComplex();
  unsigned int nD = paras.GetNDouble();
  unsigned int nI = paras.GetNInteger();

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
        paras.GetMultiComplex(0)->GetValues();
    const std::vector<std::complex<double>> v_exp =
        paras.GetMultiComplex(1)->GetValues();

    std::vector<std::complex<double>> results(v_base.size(),
                                              std::complex<double>(0., 0.));

    for (unsigned int ele = 0; ele < v_base.size(); ele++)
      results[ele] = std::pow(v_base.at(ele), v_exp.at(ele));

    out = std::shared_ptr<AbsParameter>(
        new MultiComplex(out->GetName(), results));

    break;
  } // end multi complex

  case ParType::MDOUBLE: {
    // output multi double: input must be two multi double
    if (!(nMD == 2)) {
      // TODO: exception wrong input
      return false;
    }
    // TODO: integer exponent can be calculated much faster
    const std::vector<double> v_base = paras.GetMultiDouble(0)->GetValues();
    const std::vector<double> v_exp = paras.GetMultiDouble(1)->GetValues();

    std::vector<double> results(v_base.size(), 0.);

    for (unsigned int ele = 0; ele < v_base.size(); ele++)
      results[ele] = std::pow(v_base.at(ele), v_exp.at(ele));

    out =
        std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(), results));

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
    out = std::shared_ptr<AbsParameter>(new ComplexParameter(
        out->GetName(), std::pow(tmpA->GetValue(), tmpB->GetValue())));
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
    out = std::shared_ptr<AbsParameter>(new DoubleParameter(
        out->GetName(), std::pow(tmpA->GetValue(), tmpB->GetValue())));
    break;
  } // end double

  default: {
    // TODO: exception output partype wrong
    return false;
  }

  } // end switch

  return true;
};
}
