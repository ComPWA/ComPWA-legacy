//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------

#define BOOST_TEST_MODULE Core

#include <memory>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>

#include <Core/Exceptions.hpp>
#include <Core/FunctionTree/FitParameter.hpp>
#include <Core/FunctionTree/ParameterList.hpp>
#include <Core/FunctionTree/Value.hpp>

namespace ComPWA {

using ComPWA::FunctionTree::ErrorType;
using ComPWA::FunctionTree::FitParameter;
using ComPWA::FunctionTree::ParameterList;
using ComPWA::FunctionTree::Value;

BOOST_AUTO_TEST_SUITE(ParameterTest);

BOOST_AUTO_TEST_CASE(BoundsCheck) {
  ComPWA::Logging log("", "trace");

  ComPWA::FunctionTree::FitParameter parWrong("wrongPar", 7, 1);
  BOOST_CHECK_EXCEPTION(parWrong.setBounds(20, -1), BadParameter,
                        [](const BadParameter &ex) { return true; });
  parWrong.setBounds(0, 20);
  BOOST_CHECK(parWrong.hasBounds());
}

BOOST_AUTO_TEST_CASE(SetGetCheck) {
  ComPWA::FunctionTree::FitParameter emptyInt("emptyIntPar");
  emptyInt.fixParameter(false);
  emptyInt.setValue(7);
  emptyInt.setBounds(0, 10);
  emptyInt.setError(1);
  BOOST_CHECK_CLOSE(emptyInt.value(), 7., 0.0001);
  BOOST_CHECK_CLOSE(emptyInt.bounds().first, 0., 0.0001);
  BOOST_CHECK_CLOSE(emptyInt.bounds().second, 10., 0.0001);
  BOOST_CHECK_CLOSE(emptyInt.error().first, 1., 0.0001);
}

BOOST_AUTO_TEST_CASE(FixValueCheck) {
  ComPWA::FunctionTree::FitParameter emptyFloat("emptyFloatPar");
  emptyFloat.fixParameter(false);
  emptyFloat.setValue(7.);
  emptyFloat.setBounds(0., 10.);
  emptyFloat.setError(1.);
  BOOST_CHECK(!emptyFloat.isFixed());
  emptyFloat.fixParameter(true);
  BOOST_CHECK(emptyFloat.isFixed());
  BOOST_CHECK_THROW(emptyFloat.setValue(8.), ParameterFixed);
  BOOST_CHECK_CLOSE(emptyFloat.value(), 7., 0.0001);
}

BOOST_AUTO_TEST_CASE(ConstructorCheck2) {
  Value<int> emptyInt("emptyIntPar", 1);
  ComPWA::FunctionTree::FitParameter emptyFloat("emptyFloatPar");
  ComPWA::FunctionTree::FitParameter parD("parD", 2, 0);
  ComPWA::FunctionTree::FitParameter parCopy(parD);
  ComPWA::FunctionTree::FitParameter parWrong("wrongPar", 7);
  std::shared_ptr<ComPWA::FunctionTree::FitParameter> pParInt(
      new ComPWA::FunctionTree::FitParameter("intPointerPar", 3));
  std::vector<ComPWA::FunctionTree::FitParameter> vecParInt, vecParIntCopy;
  for (unsigned int par = 0; par < 10; par++)
    vecParInt.push_back(ComPWA::FunctionTree::FitParameter(
        std::string("listPar") + std::to_string(par), par));
  vecParIntCopy = vecParInt; // copy vector

  BOOST_CHECK_EQUAL(emptyInt.value(), 1);
  BOOST_CHECK_CLOSE(emptyFloat.value(), 0., 0.0001);
  BOOST_CHECK_CLOSE(parD.value(), 2., 0.0001);
  BOOST_CHECK_CLOSE(parCopy.value(), 2., 0.0001);
  BOOST_CHECK_CLOSE(pParInt->value(), 3., 0.0001);
  for (unsigned int par = 0; par < 10; par++)
    BOOST_CHECK_CLOSE(vecParInt[par].value(), vecParIntCopy[par].value(),
                      0.0001);
}
// ----------- Testing ParameterList ------------

BOOST_AUTO_TEST_CASE(FillParameterList) {

  ParameterList list;
  for (unsigned int par = 0; par < 10; par++)
    list.addParameter(std::make_shared<ComPWA::FunctionTree::FitParameter>(
        std::string("listPar") + std::to_string(par), par, 0, 10, 1));

  std::shared_ptr<ComPWA::FunctionTree::FitParameter> dTest(
      new ComPWA::FunctionTree::FitParameter("doublePAr", 2.2));
  list.addParameter(dTest);

  auto bTest = std::make_shared<Value<int>>("IntPar", 1);
  list.addValue(bTest);

  BOOST_CHECK_EQUAL(list.numParameters(), 11);
  BOOST_CHECK_EQUAL(list.numValues(), 1);
}

BOOST_AUTO_TEST_CASE(Serialization) {
  auto shrpar = std::make_shared<ComPWA::FunctionTree::FitParameter>(
      "NNNNN", 2.5, 1.0, 3.0, 0.3);
  auto par = ComPWA::FunctionTree::FitParameter("par", 2.5, 1.0, 3.0, 0.3);
  std::ofstream ofs("paramter.xml");
  boost::archive::xml_oarchive oa(ofs, boost::archive::no_header);
  //  oa << BOOST_SERIALIZATION_NVP(*shrpar.get());
}

BOOST_AUTO_TEST_CASE(ParameterError) {
  auto par = ComPWA::FunctionTree::FitParameter("test", 1.5, 0.5);
  BOOST_CHECK_EQUAL(par.errorType(), ComPWA::ErrorType::SYM);
  BOOST_CHECK_EQUAL(par.error().first, par.error().second);
  BOOST_CHECK_EQUAL(par.hasError(), true);
  par.setError(0.2, 0.3);
  BOOST_CHECK_EQUAL(par.errorType(), ErrorType::ASYM);
  BOOST_CHECK_EQUAL(par.error().first, 0.2);
  BOOST_CHECK_EQUAL(par.error().second, 0.3);
}

BOOST_AUTO_TEST_CASE(TemplateValue) {
  auto par = std::make_shared<Value<std::vector<std::complex<double>>>>();
  ParameterList l;
  l.addValue(par);
  par->values().push_back(std::complex<double>(3, 4));
  BOOST_CHECK_EQUAL(par->value().at(0), std::complex<double>(3, 4));
}

BOOST_AUTO_TEST_SUITE_END();

} // namespace ComPWA
