// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Define Boost test module
#define BOOST_TEST_MODULE HelicityFormalism

#include <vector>

#include "Core/Logging.hpp"
#include "Core/FourMomentum.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"
#include "Physics/BuilderXML.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

#include "qft++/Vector4.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

using ComPWA::QFT::Tensor;
using ComPWA::QFT::Vector4;

BOOST_AUTO_TEST_SUITE(HelicityFormalism)

const std::string ModelConfigXML = R"####(
<ParticleList>
<Particle Name='pi0'>
	<Pid>111</Pid>
	<Parameter Class='Double' Type='Mass' Name='Mass_pi0'>
		<Value>0.1349766</Value>
		<Error>0.0000006</Error>
	</Parameter>
	<QuantumNumber Class='Spin' Type='Spin' Value='0' />
	<QuantumNumber Class='Int' Type='Charge' Value='0' />
	<QuantumNumber Class='Int' Type='Parity' Value='-1' />
	<QuantumNumber Class='Int' Type='Cparity' Value='1' />
</Particle>
<Particle Name='gamma'>
	<Pid>22</Pid>
	<Parameter Class='Double' Type='Mass' Name='mass_gamma'>
		<Value>0.</Value>
		<Fix>true</Fix>
	</Parameter>
	<QuantumNumber Class='Spin' Type='Spin' Value='1' />
	<QuantumNumber Class='Int' Type='Charge' Value='0' />
	<QuantumNumber Class='Int' Type='Parity' Value='-1' />
	<QuantumNumber Class='Int' Type='Cparity' Value='-1' />
	<QuantumNumber Class='Int' Type='Gparity' Value='1' />
</Particle>
<Particle Name='f0_980'>
	<Pid>9010221</Pid>
	<Parameter Class='Double' Type='Mass' Name='Mass_f0_980'>
		<Value>0.99</Value>
		<Fix>true</Fix>
		<Min>0.5</Min>
		<Max>1.5</Max>
		<Error>0</Error>
	</Parameter>
	<QuantumNumber Class='Spin' Type='Spin' Value='0' />
	<QuantumNumber Class='Int' Type='Charge' Value='0' />
	<QuantumNumber Class='Int' Type='Parity' Value='1' />
	<QuantumNumber Class='Int' Type='Cparity' Value='1' />
	<QuantumNumber Class='Int' Type='Gparity' Value='1' />
	<DecayInfo Type='relativisticBreitWigner'>
		<FormFactor Type='0' />
		<Parameter Class='Double' Type='Width' Name='Width_f0_980'>
			<Value>0.05</Value>
			<Fix>true</Fix>
			<Min>0.</Min>
			<Max>.5</Max>
			<Error>0</Error>
		</Parameter>
		<Parameter Class='Double' Type='MesonRadius'
			Name='Radius_f0_980'>
			<Value>1.5</Value>
			<Fix>true</Fix>
			<Min>1.0</Min>
			<Max>2.0</Max>
			<Error>0</Error>
		</Parameter>
	</DecayInfo>
</Particle>
<Particle Name='jpsi'>
	<Pid>443</Pid>
	<Parameter Class='Double' Type='Mass' Name='Mass_jpsi'>
		<Value>3.0969</Value>
		<Fix>true</Fix>
	</Parameter>
	<QuantumNumber Class='Spin' Type='Spin' Value='1' />
	<QuantumNumber Class='Int' Type='Charge' Value='0' />
	<QuantumNumber Class='Int' Type='Parity' Value='-1' />
	<QuantumNumber Class='Int' Type='Cparity' Value='-1' />
	<QuantumNumber Class='Int' Type='Gparity' Value='1' />
	<DecayInfo Type='relativisticBreitWigner'>
		<FormFactor Type='0' />
		<Parameter Class='Double' Type='Width' Name='Width_jpsi'>
			<Value>0.0000929</Value>
			<Error>0.0000028</Error>
		</Parameter>
		<Parameter Class='Double' Type='MesonRadius'
			Name='Radius_jpsi'>
			<Value>2.5</Value>
			<Fix>true</Fix>
			<Min>2.0</Min>
			<Max>3.0</Max>
		</Parameter>
	</DecayInfo>
</Particle>
<Particle Name='omega'>
	<Pid>223</Pid>
	<Parameter Class='Double' Type='Mass' Name='Mass_omega'>
		<Value>1.78265</Value>
		<Fix>true</Fix>
		<Error>0.00012</Error>
	</Parameter>
	<QuantumNumber Class='Spin' Type='Spin' Value='1' />
	<QuantumNumber Class='Int' Type='Charge' Value='0' />
	<QuantumNumber Class='Int' Type='Parity' Value='-1' />
	<QuantumNumber Class='Int' Type='Cparity' Value='-1' />
	<QuantumNumber Class='Int' Type='Gparity' Value='1' />
	<DecayInfo Type='relativisticBreitWigner'>
		<FormFactor Type='0' />
		<Parameter Class='Double' Type='Width' Name='Width_omega'>
			<Value>0.01849</Value>
			<Fix>true</Fix>
			<Min>0.0</Min>
			<Max>1.0</Max>
			<Error>0.00008</Error>
		</Parameter>
		<Parameter Class='Double' Type='MesonRadius'
			Name='Radius_omega'>
			<Value>1.5</Value>
			<Fix>true</Fix>
			<Min>1.0</Min>
			<Max>2.0</Max>
			<Error>0</Error>
		</Parameter>
	</DecayInfo>
</Particle>
</ParticleList>

<HelicityKinematics>
<PhspVolume>0.123</PhspVolume>
<InitialState>
	<Particle Name='jpsi' PositionIndex='0' />
</InitialState>
<FinalState>
	<Particle Name='pi0' Id='0' />
	<Particle Name='gamma' Id='1' />
	<Particle Name='pi0' Id='2' />
	<Particle Name='pi0' Id='3' />
</FinalState>
</HelicityKinematics>)####";

std::pair<double, double>
pawianHelicityAngles(std::vector<Vector4<double>> levels) {

  if (levels.size() < 3) {
    levels.insert(
        levels.begin(),
        Vector4<double>(std::sqrt(levels[0].Mass2() + 1.0), 0., 0., -1.0));
  }
  if (levels.size() < 4) {
    levels.insert(levels.begin(), Vector4<double>(0., 0., 0., -2.0));
  }

  Vector4<double> result = levels[3];
  Vector4<double> refTrafo = levels[1];
  Vector4<double> refRecoilTrafo = levels[0] - levels[1];
  Vector4<double> mother = levels[2];

  // boost all vectors into the mother rest frame
  result.Boost(mother);
  refTrafo.Boost(mother);
  refRecoilTrafo.Boost(mother);

  // rotate vectors so that refRecoilTrafo moves in the negative direction of
  // the z-axis
  result.RotateZ(-refTrafo.Phi());
  result.RotateY(TMath::Pi() - refTrafo.Theta());

  refRecoilTrafo.RotateZ(-refTrafo.Phi());
  refRecoilTrafo.RotateY(TMath::Pi() - refTrafo.Theta());

  // rotate around the z-axis so that refRecoil lies in the x-z plain
  result.RotateZ(M_PI - refRecoilTrafo.Phi());
  return std::make_pair(result.CosTheta(), result.Phi());
}

// ----- These functions have been copied from EvtGen and ported to qft++ -----
double evtGenDecayAngle(const Vector4<double> &p, const Vector4<double> &q,
                        const Vector4<double> &d) {
  double pd = p * d;
  double pq = p * q;
  double qd = q * d;
  double mp2 = p.Mass2();
  double mq2 = q.Mass2();
  double md2 = d.Mass2();

  double cost = (pd * mq2 - pq * qd) /
                sqrt((pq * pq - mq2 * mp2) * (qd * qd - mq2 * md2));

  return cost;
}

double mag2r3(const Vector4<double> &p0, const Vector4<double> &p1) {
  return std::pow(p0 * p1, 2) / p0.Mass2() - p1.Mass2();
}

// Calculate the 3-d dot product of 4-vectors p1 and p2 in the rest frame of
// 4-vector p0
double dotr3(const Vector4<double> &p0, const Vector4<double> &p1,
             const Vector4<double> &p2) {
  return 1 / p0.Mass2() * (p0 * p1) * (p0 * p2) - p1 * p2;
}

Tensor<double> dual(const Tensor<double> &t2) {
  Tensor<double> temp(2);
  temp(0, 0) = 0.0;
  temp(1, 1) = 0.0;
  temp(2, 2) = 0.0;
  temp(3, 3) = 0.0;

  temp(0, 1) = t2(3, 2) - t2(2, 3);
  temp(0, 2) = -t2(3, 1) + t2(1, 3);
  temp(0, 3) = t2(2, 1) - t2(1, 2);

  temp(1, 2) = -t2(3, 0) + t2(0, 3);
  temp(1, 3) = t2(2, 0) - t2(0, 2);

  temp(2, 3) = -t2(1, 0) + t2(0, 1);

  temp(1, 0) = -temp(0, 1);
  temp(2, 0) = -temp(0, 2);
  temp(3, 0) = -temp(0, 3);

  temp(2, 1) = -temp(1, 2);
  temp(3, 1) = -temp(1, 3);

  temp(3, 2) = -temp(2, 3);
  return temp;
}

// Calculate ( \vec{p1} cross \vec{p2} ) \cdot \vec{p3} in rest frame of object
double scalartripler3(const Vector4<double> &p0, const Vector4<double> &p1,
                      const Vector4<double> &p2, const Vector4<double> &p3) {
  Tensor<double> l = dual(p0 % p1) * p2;
  return -1.0 / p0.Mass() * (l * p3);
}

// Calculate phi using the given 4 vectors (all in the same frame)
double evtGenDecayAnglePhi(const Vector4<double> &z, const Vector4<double> &p,
                           const Vector4<double> &q, const Vector4<double> &d) {
  double eq = (p * q) / p.Mass();
  double ed = (p * d) / p.Mass();
  double mq = q.Mass();
  double q2 = mag2r3(p, q);
  double qd = dotr3(p, q, d);
  double zq = dotr3(p, z, q);
  double zd = dotr3(p, z, d);
  double alpha = (eq - mq) / (q2 * mq) * qd - ed / mq;

  double y = scalartripler3(p, z, q, d) + alpha * scalartripler3(p, z, q, q);
  double x = (zq * (qd + alpha * q2) - q2 * (zd + alpha * zq)) / std::sqrt(q2);

  double phi = std::atan2(y, x);

  return phi;
}

// ----- End EvtGen helper functions -----

std::pair<double, double> calculateEvtGenAngles(
    const Vector4<double> &grandparent, const Vector4<double> &parent,
    const Vector4<double> &decaying_state, const Vector4<double> &daughter) {
  return std::make_pair(
      evtGenDecayAngle(parent, decaying_state, daughter),
      evtGenDecayAnglePhi(grandparent, parent, decaying_state, daughter));
}

double calculatePhiDiff(double phi1, double phi2) {
  double phidiff = phi1 - phi2;
  if (phidiff > TMath::Pi())
    phidiff -= 2.0 * TMath::Pi();
  else if (phidiff < -TMath::Pi())
    phidiff += 2.0 * TMath::Pi();
  return phidiff;
}

/*!
 * This test case verifies the correctness of the theta and phi angles in a
 * complicated decay tree, by comparing the results with evt gen and the
 * pawian
 */
BOOST_AUTO_TEST_CASE(HelicityAnglesCorrectnessTest) {
  ComPWA::Logging log("TRACE", "output.log");

  std::stringstream ModelStringStream;
  ModelStringStream << ModelConfigXML;
  auto partL = ComPWA::readParticles(ModelStringStream);

  // Construct HelicityKinematics from XML tree
  boost::property_tree::ptree tr;
  // the boost xml parser modifies the stringstream, so that it has to be reset
  ModelStringStream.clear();
  ModelStringStream << ModelConfigXML;
  boost::property_tree::xml_parser::read_xml(ModelStringStream, tr);

  auto Kinematics = ComPWA::Physics::createHelicityKinematics(
      partL, tr.get_child("HelicityKinematics"));

  // Generate phsp sample
  ComPWA::Data::Root::RootGenerator Generator(
      Kinematics.getParticleStateTransitionKinematicsInfo());

  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(123);

  auto Sample(ComPWA::Data::generatePhsp(50, Generator, RandomGenerator));

  Vector4<double> top_vec4(0, 0, 0, 1);

  LOG(INFO) << "Loop over phsp events and comparison of angles....";
  for (auto Event : Sample.Events) {
    // convert evt to evtgen 4 vectors
    std::vector<Vector4<double>> Vectors;
    for (auto const &mom : Event.FourMomenta) {
      Vectors.push_back(Vector4<double>(mom));
    }

    Vector4<double> Level0 = Vectors[0] + Vectors[1] + Vectors[2] + Vectors[3];
    Vector4<double> Level1 = Vectors[1] + Vectors[2] + Vectors[3];
    Vector4<double> Level2 = Vectors[2] + Vectors[3];
    Vector4<double> Level3 = Vectors[2];

    // Pawian angles
    std::vector<std::pair<double, double>> PawianAngles;
    PawianAngles.push_back(pawianHelicityAngles({Level0, Level1}));
    PawianAngles.push_back(pawianHelicityAngles({Level0, Level1, Level2}));
    PawianAngles.push_back(
        pawianHelicityAngles({Level0, Level1, Level2, Level3}));

    // EvtGen angles
    std::vector<std::pair<double, double>> EvtGenAngles;
    EvtGenAngles.push_back(std::make_pair(Level1.CosTheta(), Level1.Phi()));
    EvtGenAngles.push_back(
        calculateEvtGenAngles(top_vec4, Level0, Level1, Level2));
    EvtGenAngles.push_back(
        calculateEvtGenAngles(Level0, Level1, Level2, Level3));
    Level1.Boost(Level0);
    EvtGenAngles[0] = std::make_pair(Level1.CosTheta(), Level1.Phi());

    // ComPWA angles
    using ComPWA::Physics::SubSystem;
    std::vector<SubSystem> SubSystems = {SubSystem({{1, 2, 3}, {0}}, {}, {}),
                                         SubSystem({{2, 3}, {1}}, {0}, {}),
                                         SubSystem({{2}, {3}}, {1}, {0})};

    std::vector<std::pair<double, double>> ComPwaAngles;
    for (auto SubSys : SubSystems) {
      Kinematics.registerSubSystem(SubSys);
      auto vals = Kinematics.calculateHelicityAngles(Event, SubSys);
      ComPwaAngles.push_back(std::make_pair(std::cos(vals.first), vals.second));
    }

    for (unsigned int i = 0; i < ComPwaAngles.size(); ++i) {
      BOOST_CHECK_EQUAL((float)ComPwaAngles[i].first,
                        (float)PawianAngles[i].first);
      BOOST_CHECK_EQUAL((float)ComPwaAngles[i].first,
                        (float)EvtGenAngles[i].first);
      BOOST_CHECK_EQUAL((float)ComPwaAngles[i].second,
                        (float)EvtGenAngles[i].second);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
