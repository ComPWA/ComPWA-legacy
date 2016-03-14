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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Number
#include <boost/test/unit_test.hpp>
#include <Core/Parameter.hpp>
#include <Core/ParameterList.hpp>
#include <Core/Exceptions.hpp>
#include <memory>
#include <vector>

namespace ComPWA {

BOOST_AUTO_TEST_SUITE(ParameterListSuite);

BOOST_AUTO_TEST_CASE(ConstructorCheck)
{
  std::vector<std::shared_ptr<IntegerParameter> > vecParInt;
  for(unsigned int par=0; par<10; par++)
    vecParInt.push_back(std::shared_ptr<IntegerParameter>(new  IntegerParameter(std::string("listPar"+par),par,0,10,1)));

  ParameterList testList(vecParInt);
  std::shared_ptr<DoubleParameter> dTest(new DoubleParameter("doublePAr",2.2));
  testList.AddParameter(dTest);
  std::shared_ptr<BoolParameter> bTest(new BoolParameter("boolPar",true));
  testList.AddParameter(bTest);

  BOOST_CHECK_CLOSE(testList.GetNParameter(), 12., 0.0001);
  BOOST_CHECK_CLOSE(testList.GetNDouble(), 1., 0.0001);
  BOOST_CHECK_CLOSE(testList.GetNInteger(), 10., 0.0001);
  BOOST_CHECK_CLOSE(testList.GetNBool(), 1., 0.0001);

  BOOST_CHECK_CLOSE(testList.GetIntegerParameter(4)->GetValue(), 4., 0.0001);
  BOOST_CHECK_CLOSE(testList.GetIntegerParameter(4)->GetValue(), testList.GetParameterValue(5), 0.0001);
  BOOST_CHECK(testList.GetBoolParameter(0)->GetValue());
  BOOST_CHECK_CLOSE(testList.GetBoolParameter(0)->GetValue(), testList.GetParameterValue(11), 0.0001);

  BOOST_CHECK_THROW(testList.GetBoolParameter(4), BadParameter);

  BOOST_CHECK_THROW(testList.SetParameterValue(2,5.), BadParameter);
  testList.SetParameterValue(0,1.1);
  BOOST_CHECK_CLOSE(testList.GetDoubleParameter(0)->GetValue(), 1.1, 0.0001);
}

BOOST_AUTO_TEST_SUITE_END();

} /* namespace ComPWA */
