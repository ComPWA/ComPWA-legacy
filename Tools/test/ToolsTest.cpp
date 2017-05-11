
#define BOOST_TEST_MODULE ToolsTest

#include <iostream>
#include <boost/test/unit_test.hpp>
#include "Tools/Integration.hpp"

using namespace ComPWA::Tools;

BOOST_AUTO_TEST_SUITE(ToolsTest)

BOOST_AUTO_TEST_CASE(IntegrationTest){
  auto gauss = testGauss();
  auto intAlg = IntegralByQuadrature<testGauss>(gauss, std::pair<double,double>(-10,10));
  std::cout<<intAlg.Integral(1000)<<std::endl;
  BOOST_CHECK( std::abs(intAlg.Integral() - 1.) < 1e-3);
  
  std::cout<<gauss(0)<<std::endl;
};
BOOST_AUTO_TEST_SUITE_END()
