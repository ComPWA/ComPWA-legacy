#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Number
#include <boost/test/unit_test.hpp>
#include <PWAGenericPar.hpp>
#include <memory>
#include <vector>

BOOST_AUTO_TEST_SUITE(GenericParSuite);

BOOST_AUTO_TEST_CASE(BoundsCheck)
{
  PWAGenericPar<int> parWrong(7,10,5,1);
  BOOST_CHECK(!parWrong.HasBounds());
  parWrong.SetTMaxValue(-1);
  BOOST_CHECK(!parWrong.HasBounds());
  parWrong.SetTMaxValue(6);
  BOOST_CHECK(!parWrong.HasBounds());
  parWrong.SetTMinMax(10,5);
  BOOST_CHECK(!parWrong.HasBounds());
  parWrong.SetTMinMax(5,10);
  BOOST_CHECK(parWrong.HasBounds());
}

BOOST_AUTO_TEST_CASE(SetGetCheck)
{
  PWAGenericPar<int> emptyInt;
  emptyInt.SetValue(7); emptyInt.SetMinMax(0,10); emptyInt.SetError(1);
  BOOST_CHECK_CLOSE(emptyInt.GetValue(), 7., 0.0001);
  BOOST_CHECK_CLOSE(emptyInt.GetMinValue(), 0., 0.0001);
  BOOST_CHECK_CLOSE(emptyInt.GetMaxValue(), 10., 0.0001);
  BOOST_CHECK_CLOSE(emptyInt.GetError(), 1., 0.0001);
}

BOOST_AUTO_TEST_CASE(ConstructorCheck)
{
  PWAGenericPar<int> emptyInt;
  PWAGenericPar<double> emptyFloat;
  PWAGenericPar<int> parInt(2,0,5,1);
  PWAGenericPar<int> parCopy(parInt);
  PWAGenericPar<int> parWrong(7,10,5,1);
  std::shared_ptr<PWAGenericPar<int> > pParInt(new PWAGenericPar<int>(3,0,5,1));
  std::vector<PWAGenericPar<int> > vecParInt, vecParIntCopy;
  for(unsigned int par=0; par<10; par++)
    vecParInt.push_back(PWAGenericPar<int>(par,0,10,1));
  vecParIntCopy = vecParInt; //copy vector

  BOOST_CHECK_CLOSE(emptyInt.GetValue(), 0., 0.0001);
  BOOST_CHECK_CLOSE(emptyFloat.GetValue(), 0., 0.0001);
  BOOST_CHECK_CLOSE(parInt.GetValue(), 2., 0.0001);
  BOOST_CHECK_CLOSE(parCopy.GetValue(), 2., 0.0001);
  BOOST_CHECK(!parWrong.HasBounds());
  BOOST_CHECK_CLOSE(pParInt->GetValue(), 3., 0.0001);
  for(unsigned int par=0; par<10; par++)
    BOOST_CHECK_CLOSE(vecParInt[par].GetValue(), vecParIntCopy[par].GetValue(), 0.0001);

}

BOOST_AUTO_TEST_SUITE_END();
