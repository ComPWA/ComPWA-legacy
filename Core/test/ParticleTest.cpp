
#define BOOST_TEST_MODULE Core

#include <memory>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <Core/Exceptions.hpp>
#include <Core/Particle.hpp>

namespace ComPWA {

BOOST_AUTO_TEST_SUITE(ParticleTest);

BOOST_AUTO_TEST_CASE(FourMomentum) {
  ComPWA::FourMomentum p4(1,2,3,4);
  ComPWA::FourMomentum p4B(1,2,3,4);
  BOOST_CHECK_EQUAL(p4, p4B);
  BOOST_CHECK_EQUAL(p4.invMassSq(), 2.0);
  BOOST_CHECK_EQUAL(ComPWA::FourMomentum::threeMomentumSq(p4), 14.0);
  
  p4B.setValue( std::array<double,4>{{1,2,3,5}} );
  auto pTot = p4B+p4; //(2,4,6,9)
  BOOST_CHECK_EQUAL(pTot.invMassSq(), 25.0);
  
  ComPWA::FourMomentum pSum;
  pSum += p4;
  pSum += p4B;
  BOOST_CHECK_EQUAL(pSum.invMassSq(), 25.0);
}
  
BOOST_AUTO_TEST_CASE(Particle) {
  ComPWA::Particle part(1,2,3,4);
  BOOST_CHECK_EQUAL(part.massSq(), 2.0);
}

BOOST_AUTO_TEST_SUITE_END();

} /* namespace ComPWA */
