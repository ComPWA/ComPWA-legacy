// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Create a global string which contains kind of PDG database in ComPWA format
///

///
/// Gobal particle list. Incomplete list based on PDG2016.
/// The string can be converted to a PartList whereever needed:
/// \code
///   #include "Physics/ParticleList.hpp"
///   # convert to stream
///   std::stringstream str;
///   str << defaultParticleList;
///   # fill property tree
///   boost::property_tree::ptree tr;
///   boost::property_tree::xml_parser::read_xml(str, tr);
///   # fill PartList
///   auto partL = std::make_shared<ComPWA::PartList>();
///   ReadParticles(partL, tr);
/// \endcode
/// We use C++11 raw string literal here.
///
namespace ComPWA {
namespace Physics {

const std::string defaultParticleList = R"####(
<ParticleList>
	<Particle Name="gamma">
		<Pid>22</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_gamma">
			<Value>0.0</Value>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="1"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="1"/>
	</Particle>

	# Leptons
	<Particle Name="e+">
		<Pid>11</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_electron">
			<Value>0.0005109989461</Value>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0.5"/>
		<QuantumNumber Class="Int" Type="Charge" Value="1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
	</Particle>
	<Particle Name="e-">
		<Pid>-11</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_electron">
			<Value>0.0005109989461</Value>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0.5"/>
		<QuantumNumber Class="Int" Type="Charge" Value="-1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="+1"/>
	</Particle>
	<Particle Name="mu+">
		<Pid>13</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_muon">
			<Value>0.1056583745</Value>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0.5"/>
		<QuantumNumber Class="Int" Type="Charge" Value="1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
	</Particle>
	<Particle Name="mu-">
		<Pid>-13</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_muon">
			<Value>0.1056583745</Value>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0.5"/>
		<QuantumNumber Class="Int" Type="Charge" Value="-1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="+1"/>
	</Particle>
	<Particle Name="tau+">
		<Pid>15</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_tau">
			<Value>1.77686</Value>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0.5"/>
		<QuantumNumber Class="Int" Type="Charge" Value="1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
	</Particle>
	<Particle Name="tau-">
		<Pid>-15</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_tau">
			<Value>1.77686</Value>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0.5"/>
		<QuantumNumber Class="Int" Type="Charge" Value="-1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="+1"/>
	</Particle>

	# Light flavoured mesons
	<Particle Name="pi0">
		<Pid>111</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_neutralPion">
			<Value>0.1349766</Value>
			<Error>0.000006</Error>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="1" Projection="0"/>
	</Particle>
	<Particle Name="pi+">
		<Pid>211</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_chargedPion">
			<Value>0.13957018</Value>
			<Error>0.00000035</Error>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="+1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="1" Projection="1"/>
	</Particle>
	<Particle Name="pi-">
		<Pid>-211</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_chargedPion">
			<Value>0.13957018</Value>
			<Error>0.00000035</Error>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="-1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="1" Projection="-1"/>
	</Particle>
	<Particle Name="eta">
		<Pid>221</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_eta">
			<Value>0.547862</Value>
			<Error>0.000017</Error>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="1"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Width" Name="Width_eta">
				<Value>0.00000131</Value>
				<Error>0.00000005</Error>
				<Fix>true</Fix>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_eta">
				<Value>2.5</Value>
				<Fix>true</Fix>
				<Min>2.0</Min>
				<Max>3.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="rho(770)0">
		<Pid>113</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_rho770">
			<Value>0.7526</Value>
			<Error>0.00025</Error>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="1" Projection="0"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Width" Name="Width_rho">
				<Value>0.000001491</Value>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_rho">
				<Value>2.5</Value>
				<Fix>true</Fix>
				<Min>2.0</Min>
				<Max>3.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="rho(770)+">
		<Pid>213</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_rho770">
			<Value>0.7526</Value>
			<Error>0.00025</Error>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="+1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="1" Projection="+1"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Width" Name="Width_rho">
				<Value>0.000001491</Value>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_rho">
				<Value>2.5</Value>
				<Fix>true</Fix>
				<Min>2.0</Min>
				<Max>3.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="rho(770)-">
		<Pid>-213</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_rho770">
			<Value>0.7526</Value>
			<Error>0.00025</Error>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="-1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="1" Projection="-1"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Width" Name="Width_rho">
				<Value>0.000001491</Value>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_rho">
				<Value>2.5</Value>
				<Fix>true</Fix>
				<Min>2.0</Min>
				<Max>3.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="omega(782)">
		<Pid>223</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_omega782">
			<Value>0.78265</Value>
			<Error>0.00012</Error>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="1"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="0"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Width" Name="Width_omega782">
				<Value>0.000001491</Value>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_omega782">
				<Value>2.5</Value>
				<Fix>true</Fix>
				<Min>2.0</Min>
				<Max>3.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="a0(980)0">
		<Pid>9000111</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_a0(980)">
			<Value>0.994</Value>
			<Error>0.001</Error>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="1"/>
		<DecayInfo Type="flatte">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Coupling" Name="gKK_a0(980)">
				<Value>3.121343843602647</Value>
				<Error>0.001</Error>
				<Fix>false</Fix>
				<ParticleA>K+</ParticleA>
				<ParticleB>K-</ParticleB>
			</Parameter>
			<Parameter Class="Double" Type="Coupling" Name="gEtaPi_a0(980)">
				<Value>2.66</Value>
				<Error>0.001</Error>
				<Fix>true</Fix>
				<ParticleA>eta</ParticleA>
				<ParticleB>pi0</ParticleB>
			</Parameter>
			<Parameter Class="Double" Type="Coupling" Name="gKK_a0(980)">
				<Value>3.121343843602647</Value>
				<Error>0.001</Error>
				<Fix>false</Fix>
				<ParticleA>K_S0</ParticleA>
				<ParticleB>K_S0</ParticleB>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_a0(980)">
				<Value>1.5</Value>
				<Fix>true</Fix>
				<Min>1.0</Min>
				<Max>2.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="a0(980)+">
		<Pid>9000211</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_a0(980)">
			<Value>0.994</Value>
			<Error>0.001</Error>
			<Fix>true</Fix>
		</Parameter>
		<Charge>0</Charge>
		<Spin>0</Spin>
		<Parity>+1</Parity>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="+1"/>
		<DecayInfo Type="flatte">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Coupling" Name="gKK_a0(980)">
				<Value>3.121343843602647</Value>
				<Error>0.001</Error>
				<Fix>false</Fix>
				<ParticleA>K_S0</ParticleA>
				<ParticleB>K+</ParticleB>
			</Parameter>
			<Parameter Class="Double" Type="Coupling" Name="gEtaPi_a0(980)">
				<Value>2.66</Value>
				<Error>0.001</Error>
				<Fix>true</Fix>
				<ParticleA>eta</ParticleA>
				<ParticleB>pi0</ParticleB>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_a0(980)">
				<Value>1.5</Value>
				<Fix>true</Fix>
				<Min>1.0</Min>
				<Max>2.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="a0(980)-">
		<Pid>9000211</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_a0(980)">
			<Value>0.994</Value>
			<Error>0.001</Error>
			<Fix>true</Fix>
		</Parameter>
		<Charge>0</Charge>
		<Spin>0</Spin>
		<Parity>+1</Parity>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="+1"/>
		<DecayInfo Type="flatte">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Coupling" Name="gKK_a0(980)">
				<Value>3.121343843602647</Value>
				<Error>0.001</Error>
				<Fix>false</Fix>
				<ParticleA>K_S0</ParticleA>
				<ParticleB>K+</ParticleB>
			</Parameter>
			<Parameter Class="Double" Type="Coupling" Name="gEtaPi_a0(980)">
				<Value>2.66</Value>
				<Error>0.001</Error>
				<Fix>true</Fix>
				<ParticleA>eta</ParticleA>
				<ParticleB>pi0</ParticleB>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_a0(980)">
				<Value>1.5</Value>
				<Fix>true</Fix>
				<Min>1.0</Min>
				<Max>2.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="phi(1020)">
		<Pid>333</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_phi(1020)">
			<Value>1.019461</Value>
			<Error>0.000019</Error>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="1"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="1" />
			<Parameter Class="Double" Type="Width" Name="Width_phi(1020)">
				<Value>0.004266</Value>
				<Error>0.000031</Error>
				<Fix>true</Fix>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_phi(1020)">
				<Value>1.5</Value>
				<Fix>true</Fix>
				<Min>1.0</Min>
				<Max>2.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>

	# Strange mesons
	<Particle Name="K-">
		<Pid>-321</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_chargedKaon">
			<Value>0.493677</Value>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="-1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="0.5" Projection="-0.5"/>
	</Particle>
	<Particle Name="K+">
		<Pid>321</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_chargedKaon">
			<Value>0.493677</Value>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="0.5" Projection="0.5"/>
	</Particle>
	<Particle Name="K_S0">
		<Pid>310</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_neutralKaon">
			<Value>0.497614</Value>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="0.5" Projection="0.5"/>
	</Particle>
	<Particle Name="K_L0">
		<Pid>130</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_neutralKaon">
			<Value>0.497614</Value>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="0.5" Projection="-0.5"/>
	</Particle>

	# Heavy flavour mesons
	<Particle Name="D+">
		<Pid>411</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_chargedD">
			<Value>1.86958</Value>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="+1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="0.5"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Width" Name="Width_chargedD">
				<Value>6.33E-13</Value>
				<Fix>true</Fix>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_chargedD">
				<Value>2.5</Value>
				<Fix>true</Fix>
				<Min>2.0</Min>
				<Max>3.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="D-">
		<Pid>-411</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_chargedD">
			<Value>1.86958</Value>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="-1"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="0.5"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Width" Name="Width_chargedD">
				<Value>6.33E-13</Value>
				<Fix>true</Fix>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_chargedD">
				<Value>2.5</Value>
				<Fix>true</Fix>
				<Min>2.0</Min>
				<Max>3.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="D0">
		<Pid>421</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_neutralD">
			<Value>1.86483</Value>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="0"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Spin" Type="IsoSpin" Value="0.5"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Width" Name="Width_neutralD">
				<Value>1.605E-12</Value>
				<Fix>true</Fix>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_neutralD">
				<Value>2.5</Value>
				<Fix>true</Fix>
				<Min>2.0</Min>
				<Max>3.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="J/psi">
		<Pid>443</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_jpsi">
			<Value>3.096900</Value>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="1"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="-1"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Class="Double" Type="Width" Name="Width_jpsi">
				<Value>9.29E-05</Value>
				<Fix>true</Fix>
			</Parameter>
			<Parameter Class="Double" Type="MesonRadius" Name="Radius_jpsi">
				<Value>2.5</Value>
				<Fix>true</Fix>
				<Min>2.0</Min>
				<Max>3.0</Max>
			</Parameter>
		</DecayInfo>
	</Particle>

</ParticleList>
)####";
}
}
