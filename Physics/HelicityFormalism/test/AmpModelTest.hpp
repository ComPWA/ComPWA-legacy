//
//  AmpModelTest.hpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 04.12.17.
//

#ifndef AmpModelTest_h
#define AmpModelTest_h

const std::string HelicityTestParticles = R"####(
<ParticleList>
  <Particle Name='K+'>
    <Pid>321</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_chargedKaon'>
      <Value>0.493677</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='1'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Spin' Type='IsoSpin' Value='0.5' Projection='0.5'/>
  </Particle>
  <Particle Name='K-'>
    <Pid>-321</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_chargedKaon'>
      <Value>0.493677</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='-1'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Spin' Type='IsoSpin' Value='0.5' Projection='-0.5'/>
  </Particle>
  <Particle Name='pi0'>
    <Pid>111</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_pi0'>
      <Value>0.1349766</Value>
      <Error>0.0000006</Error>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='1'/>
  </Particle>
  <Particle Name='eta'>
    <Pid>221</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_eta'>
      <Value>0.547862</Value>
      <Error>0.000018</Error>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='1'/>
  </Particle>
  <Particle Name='K_S0'>
    <Pid>310</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_neutralKaon'>
      <Value>0.497614</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Spin' Type='IsoSpin' Value='0.5' Projection='0.5'/>
  </Particle>
  <Particle Name='D0'>
    <Pid>421</Pid>
    <Parameter Type='Mass' Name='Mass_D0'>
      <Value>1.86484</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='0' />
      <Parameter Class='Double' Type='Mass' Name='Width_D0'>
        <Value>0.000623</Value>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_D0'>
        <Value>2.5</Value>
        <Fix>true</Fix>
        <Min>2.0</Min>
        <Max>3.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name='phi(1020)'>
    <Pid>333</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_phi(1020)'>
      <Value>1.019461</Value>
      <Error>0.000019</Error>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='1' />
      <Parameter Class='Double' Type='Width' Name='Width_phi(1020)'>
        <Value>0.004266</Value>
        <Error>0.000031</Error>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_phi(1020)'>
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name='a0(980)0'>
    <Pid>9000111</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_a0(980)'>
      <Value>0.994</Value>
      <Error>0.001</Error>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='1'/>
    <DecayInfo Type='flatte'>
      <FormFactor Type='0' />
      <Parameter Class='Double' Type='Coupling' Name='gKK_a0(980)'>
        <Value>3.121343843602647</Value>
        <Error>0.001</Error>
        <Fix>false</Fix>
        <ParticleA>K+</ParticleA>
        <ParticleB>K-</ParticleB>
      </Parameter>
      <Parameter Class='Double' Type='Coupling' Name='gEtaPi_a0(980)'>
        <Value>2.66</Value>
        <Error>0.001</Error>
        <Fix>true</Fix>
        <ParticleA>eta</ParticleA>
        <ParticleB>pi0</ParticleB>
      </Parameter>
      <Parameter Class='Double' Type='Coupling' Name='gKK_a0(980)'>
        <Value>3.121343843602647</Value>
        <Error>0.001</Error>
        <Fix>false</Fix>
        <ParticleA>K_S0</ParticleA>
        <ParticleB>K_S0</ParticleB>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_a0(980)'>
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name='a0(980)+'>
    <Pid>9000211</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_a0(980)'>
      <Value>0.994</Value>
      <Error>0.001</Error>
      <Fix>true</Fix>
    </Parameter>
    <Charge>0</Charge>
    <Spin>0</Spin>
    <Parity>+1</Parity>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='+1'/>
    <DecayInfo Type='flatte'>
      <FormFactor Type='0' />
      <Parameter Class='Double' Type='Coupling' Name='gKK_a0(980)'>
        <Value>3.121343843602647</Value>
        <Error>0.001</Error>
        <Fix>false</Fix>
        <ParticleA>K_S0</ParticleA>
        <ParticleB>K+</ParticleB>
      </Parameter>
      <Parameter Class='Double' Type='Coupling' Name='gEtaPi_a0(980)'>
        <Value>2.66</Value>
        <Error>0.001</Error>
        <Fix>true</Fix>
        <ParticleA>eta</ParticleA>
        <ParticleB>pi0</ParticleB>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_a0(980)'>
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name='Bkgphi(1020)'>
    <Pid>3339999</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_phi(1020)'>
      <Value>1.019461</Value>
      <Error>0.000019</Error>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='1' />
      <Parameter Class='Double' Type='Width' Name='Width_phi(1020)'>
        <Value>0.004266</Value>
        <Error>0.000031</Error>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_phi(1020)'>
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name='something'>
    <Pid>123456</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_something'>
      <Value>1.024654338170585</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='1' />
      <Parameter Class='Double' Type='Width' Name='Width_something'>
        <Value>0.01909658663476588</Value>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_something'>
        <Value>1.5</Value>
        <Fix>true</Fix>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name='a2(1320)-'>
    <Pid>215</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_a2(1320)'>
      <Value>1.3181</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='2'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='1' />
      <Parameter Class='Double' Type='Width' Name='Width_a2(1320)'>
        <Value>0.1098</Value>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_a2(1320)'>
        <Value>1.5</Value>
        <Fix>true</Fix>
      </Parameter>
    </DecayInfo>
  </Particle>

  <Particle Name='gamma'>
    <Pid>22</Pid>
    <Parameter Class='Double' Type='Mass' Name='mass_gamma'>
      <Value>0.</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
  </Particle>
  <Particle Name='jpsi'>
    <Pid>443</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_jpsi'>
      <Value>3.0969</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='0' />
      <Parameter Class='Double' Type='Width' Name='Width_jpsi'>
        <Value>0.0000929</Value>
        <Error>0.0000028</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_jpsi'>
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
      <Value>0.78265</Value>
      <Fix>true</Fix>
      <Min>0.5</Min>
      <Max>1.5</Max>
      <Error>0.00012</Error>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='0' />
      <Parameter Class='Double' Type='Width' Name='Width_omega'>
        <Value>0.01849</Value>
        <Fix>true</Fix>
        <Min>0.0</Min>
        <Max>1.0</Max>
        <Error>0.00008</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_omega'>
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
        <Error>0</Error>
      </Parameter>
    </DecayInfo>
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
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='0' />
      <Parameter Class='Double' Type='Width' Name='Width_f0_980'>
        <Value>0.05</Value>
        <Fix>true</Fix>
        <Min>0.</Min>
        <Max>.5</Max>
        <Error>0</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_f0_980'>
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
        <Error>0</Error>
      </Parameter>
    </DecayInfo>
  </Particle>
</ParticleList>
)####";

const std::string HelicityTestKinematics = R"####(
<HelicityKinematics>
  <PhspVolume>0.123</PhspVolume>
  <InitialState>
    <Particle Name='jpsi' PositionIndex='0'/>
  </InitialState>
  <FinalState>
    <Particle Name='pi0' Id='0'/>
    <Particle Name='gamma' Id='1'/>
    <Particle Name='pi0' Id='2'/>
  </FinalState>
</HelicityKinematics>
)####";

const std::string HelicityTestModel = R"####(
<Intensity Class='Incoherent' Name='jpsiToPi0Pi0Gamma_inc'>
  <Parameter Class='Double' Type='Strength' Name="Strength_jpsiToPi0Pi0Gamm_inc">
    <Value>1</Value>
    <Fix>true</Fix>
  </Parameter>
  <Intensity Class='Coherent' Name='jpsiToPi0Pi0Gamma'>
    <Parameter Class='Double' Type='Strength' Name='Strength_jpsiToPi0Pi0Gamma'>
      <Value>0.99</Value>
      <Fix>true</Fix>
    </Parameter>
    <Amplitude Class='SequentialPartialAmplitude' Name='omega'>
      <Parameter Class='Double' Type='Magnitude' Name='Magnitude_sosososos'>
        <Value>1.0</Value>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type='Phase' Name='Phase_sosososos'>
        <Value>0.0</Value>
        <Fix>true</Fix>
      </Parameter>
      <PartialAmplitude Class="HelicityDecay" Name='jpsitoOmegaPi0'>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsitoOmegaPi0'>
          <Value>1.0</Value>
          <Fix>true</Fix>
          <Min>0.5</Min>
          <Max>1.5</Max>
          <Error>0</Error>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsitoOmegaPi0'>
          <Value>1.0</Value>
          <Fix>true</Fix>
          <Min>0.5</Min>
          <Max>1.5</Max>
          <Error>0</Error>
        </Parameter>
        <DecayParticle Name='jpsi' Helicity='+1' />
        <DecayProducts>
          <Particle Name='omega' FinalState='0 1' Helicity='+1' />
          <Particle Name='pi0' FinalState='2' Helicity='0' />
        </DecayProducts>
      </PartialAmplitude>
      <PartialAmplitude Class="HelicityDecay" Name="omegatoPi0G">
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_omegaToPi0Gamma'>
          <Value>1.0</Value>
          <Fix>true</Fix>
          <Min>0.5</Min>
          <Max>1.5</Max>
          <Error>0</Error>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_omegaToPi0Gamma'>
          <Value>1.0</Value>
          <Fix>true</Fix>
          <Min>0.5</Min>
          <Max>1.5</Max>
          <Error>0</Error>
        </Parameter>
        <DecayParticle Name='omega' Helicity='+1' />
        <RecoilSystem FinalState='2' />
        <DecayProducts>
          <Particle Name='gamma' FinalState='1' Helicity='+1' />
          <Particle Name='pi0' FinalState='0' Helicity='0' />
        </DecayProducts>
      </PartialAmplitude>
    </Amplitude>
  </Intensity>
</Intensity>
)####";

const std::string SeqPartialAmplitudeTestModel = R"####(
    <Amplitude Class='SequentialPartialAmplitude' Name='omega'>
      <Parameter Class='Double' Type='Magnitude' Name='Magnitude_sosososos'>
        <Value>1.0</Value>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type='Phase' Name='Phase_sosososos'>
        <Value>0.0</Value>
        <Fix>true</Fix>
      </Parameter>
      <PartialAmplitude Class="HelicityDecay" Name='jpsitoOmegaPi0'>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsitoOmegaPi0'>
          <Value>1.0</Value>
          <Fix>true</Fix>
          <Min>0.5</Min>
          <Max>1.5</Max>
          <Error>0</Error>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsitoOmegaPi0'>
          <Value>1.0</Value>
          <Fix>true</Fix>
          <Min>0.5</Min>
          <Max>1.5</Max>
          <Error>0</Error>
        </Parameter>
        <DecayParticle Name='jpsi' Helicity='+1' />
        <DecayProducts>
          <Particle Name='omega' FinalState='0 1' Helicity='+1' />
          <Particle Name='pi0' FinalState='2' Helicity='0' />
        </DecayProducts>
      </PartialAmplitude>
      <PartialAmplitude Class="HelicityDecay" Name="omegatoPi0G">
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_omegaToPi0Gamma'>
          <Value>1.0</Value>
          <Fix>true</Fix>
          <Min>0.5</Min>
          <Max>1.5</Max>
          <Error>0</Error>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_omegaToPi0Gamma'>
          <Value>1.0</Value>
          <Fix>true</Fix>
          <Min>0.5</Min>
          <Max>1.5</Max>
          <Error>0</Error>
        </Parameter>
        <DecayParticle Name='omega' Helicity='+1' />
        <RecoilSystem FinalState='2' />
        <DecayProducts>
          <Particle Name='gamma' FinalState='1' Helicity='+1' />
          <Particle Name='pi0' FinalState='0' Helicity='0' />
        </DecayProducts>
      </PartialAmplitude>
    </Amplitude>
)####";

const std::string PartialAmplitudeTestModel = R"####(
      <PartialAmplitude Class="HelicityDecay" Name="omegatoPi0G">
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_omegaToPi0Gamma'>
          <Value>1.0</Value>
          <Fix>true</Fix>
          <Min>0.5</Min>
          <Max>1.5</Max>
          <Error>0</Error>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_omegaToPi0Gamma'>
          <Value>1.0</Value>
          <Fix>true</Fix>
          <Min>0.5</Min>
          <Max>1.5</Max>
          <Error>0</Error>
        </Parameter>
        <DecayParticle Name='omega' Helicity='+1' />
        <RecoilSystem FinalState='2' />
        <DecayProducts>
          <Particle Name='gamma' FinalState='1' Helicity='+1' />
          <Particle Name='pi0' FinalState='0' Helicity='0' />
        </DecayProducts>
      </PartialAmplitude>
)####";

#endif
