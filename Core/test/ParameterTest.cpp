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

#include <boost/test/unit_test.hpp>
#include <Core/Parameter.hpp>
#include <Core/Exceptions.hpp>
#include <Core/ParameterList.hpp>

namespace ComPWA {

BOOST_AUTO_TEST_SUITE(ParameterTest);

BOOST_AUTO_TEST_CASE(BoundsCheck) {
  IntegerParameter parWrong("wrongPar", 7, 10, 5, 1);
  BOOST_CHECK(!parWrong.HasBounds());
  parWrong.SetMaxValue(-1);
  BOOST_CHECK(!parWrong.HasBounds());
  parWrong.SetMaxValue(6);
  BOOST_CHECK(!parWrong.HasBounds());
  parWrong.SetMinMax(10, 5);
  BOOST_CHECK(!parWrong.HasBounds());
  parWrong.SetMinMax(5, 10);
  BOOST_CHECK(parWrong.HasBounds());
}

BOOST_AUTO_TEST_CASE(SetGetCheck) {
  IntegerParameter emptyInt("emptyIntPar");
  emptyInt.SetValue(7);
  emptyInt.SetMinMax(0, 10);
  emptyInt.SetError(1);
  BOOST_CHECK_CLOSE(emptyInt.GetValue(), 7., 0.0001);
  BOOST_CHECK_CLOSE(emptyInt.GetMinValue(), 0., 0.0001);
  BOOST_CHECK_CLOSE(emptyInt.GetMaxValue(), 10., 0.0001);
  BOOST_CHECK_CLOSE(emptyInt.GetError(), 1., 0.0001);
}

BOOST_AUTO_TEST_CASE(FixValueCheck) {
  DoubleParameter emptyFloat("emptFloatPar");
  emptyFloat.SetValue(7.);
  emptyFloat.SetMinMax(0., 10.);
  emptyFloat.SetError(1.);
  BOOST_CHECK(!emptyFloat.IsFixed());
  emptyFloat.SetParameterFixed();
  BOOST_CHECK(emptyFloat.IsFixed());
  BOOST_CHECK_THROW(emptyFloat.SetValue(8.), ParameterFixed);
  BOOST_CHECK_CLOSE(emptyFloat.GetValue(), 7., 0.0001);
  BOOST_CHECK_EQUAL(1, 2);
}

BOOST_AUTO_TEST_CASE(ConstructorCheck2) {
  IntegerParameter emptyInt("emptyIntPar");
  DoubleParameter emptyFloat("emptyFloatPar");
  IntegerParameter parInt("intPar", 2, 0, 5, 1);
  IntegerParameter parCopy(parInt);
  IntegerParameter parWrong("wrongPar", 7, 10, 5, 1);
  std::shared_ptr<IntegerParameter> pParInt(
      new IntegerParameter("intPointerPar", 3, 0, 5, 1));
  std::vector<DoubleParameter> vecParInt, vecParIntCopy;
  for (unsigned int par = 0; par < 10; par++)
    vecParInt.push_back(DoubleParameter(
        std::string("listPar") + std::to_string(par), par, 0, 10, 1));
  vecParIntCopy = vecParInt; // copy vector

  BOOST_CHECK_CLOSE(emptyInt.GetValue(), 0., 0.0001);
  BOOST_CHECK_CLOSE(emptyFloat.GetValue(), 0., 0.0001);
  BOOST_CHECK_CLOSE(parInt.GetValue(), 2., 0.0001);
  BOOST_CHECK_CLOSE(parCopy.GetValue(), 2., 0.0001);
  BOOST_CHECK(!parWrong.HasBounds());
  BOOST_CHECK_CLOSE(pParInt->GetValue(), 3., 0.0001);
  for (unsigned int par = 0; par < 10; par++)
    BOOST_CHECK_CLOSE(vecParInt[par].GetValue(), vecParIntCopy[par].GetValue(),
                      0.0001);
}
  /* ----------- Testing ParameterList ------------ */

BOOST_AUTO_TEST_CASE(ConstructorCheck) {
  std::vector<std::shared_ptr<IntegerParameter>> vecParInt;
  for (unsigned int par = 0; par < 10; par++)
    vecParInt.push_back(std::make_shared<IntegerParameter>(
        std::string("listPar") + std::to_string(par), par, 0, 10, 1));

  ParameterList testList(vecParInt);
  std::shared_ptr<DoubleParameter> dTest(new DoubleParameter("doublePAr", 2.2));
  testList.AddParameter(dTest);
  std::shared_ptr<BoolParameter> bTest(new BoolParameter("boolPar", true));
  testList.AddParameter(bTest);

  BOOST_CHECK_CLOSE(testList.GetNParameter(), 12., 0.0001);
  BOOST_CHECK_CLOSE(testList.GetNDouble(), 1., 0.0001);
  BOOST_CHECK_CLOSE(testList.GetNInteger(), 10., 0.0001);
  BOOST_CHECK_CLOSE(testList.GetNBool(), 1., 0.0001);

  BOOST_CHECK_CLOSE(testList.GetIntegerParameter(4)->GetValue(), 4., 0.0001);
  BOOST_CHECK(testList.GetBoolParameter(0)->GetValue());

  BOOST_CHECK_THROW(testList.GetBoolParameter(4), BadParameter);

  BOOST_CHECK_THROW(testList.SetParameterValue(2, 5.), BadParameter);
  testList.SetParameterValue(0, 1.1);
  BOOST_CHECK_CLOSE(testList.GetDoubleParameter(0)->GetValue(), 1.1, 0.0001);
}

BOOST_AUTO_TEST_SUITE_END();

} /* namespace ComPWA */
