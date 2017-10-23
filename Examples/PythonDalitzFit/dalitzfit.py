#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Dalitz_ext import *

import numpy as np



############################# Configuration ##############

amplitudeModel = '''
<IncoherentIntensity Name="jpsiGammaPiPi_inc">
	<CoherentIntensity Name="jpsiGammaPiPi">
  	<Amplitude Name="f2(1270)">
			<Parameter Type="Magnitude"	Name="Magnitude_f2">
				<Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>2.0</Max>
				<Fix>false</Fix>
			</Parameter>
			<Parameter Type="Phase" Name="Phase_f2">
				<Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
				<Fix>false</Fix>
			</Parameter>
			<Resonance Name="f2ToPiPi">
				<DecayParticle Name="f2(1270)" Helicity="0"/>
				<SubSystem>
					<RecoilSystem FinalState="0" />
					<DecayProducts>
						<Particle Name="pi0" FinalState="1"  Helicity="0"/>
						<Particle Name="pi0" FinalState="2"  Helicity="0"/>
					</DecayProducts>
				</SubSystem>
			</Resonance>
		</Amplitude>
		<Amplitude Name="myAmp">
			<Parameter Type="Magnitude"	Name="Magnitude_my">
				<Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>2.0</Max>
				<Fix>true</Fix>
			</Parameter>
			<Parameter Type="Phase" Name="Phase_my`">
				<Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
				<Fix>true</Fix>
			</Parameter>
			<Resonance Name="MyResToPiPi">
				<DecayParticle Name="myRes" Helicity="0"/>
				<SubSystem>
					<RecoilSystem FinalState="0" />
					<DecayProducts>
						<Particle Name="pi0" FinalState="1"  Helicity="0"/>
						<Particle Name="pi0" FinalState="2"  Helicity="0"/>
					</DecayProducts>
				</SubSystem>
			</Resonance>
		</Amplitude>
	</CoherentIntensity>
</IncoherentIntensity>
'''

myParticles = '''
<ParticleList>
	<Particle Name="f2(1270)">
		<Pid>225</Pid>
		<Parameter Type="Mass" Name="Mass_f2(1270)">
			<Value>1.2755</Value>
			<Error>8.0E-04</Error>
      <Min>0.1</Min>
      <Max>2.0</Max>
      <Fix>false</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="2"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="+1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Type="Width" Name="Width_f2(1270)">
				<Value>0.1867</Value>
			</Parameter>
			<Parameter Type="MesonRadius" Name="Radius_rho">
				<Value>2.5</Value>
				<Fix>true</Fix>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="myRes">
		<Pid>999999</Pid>
		<Parameter Type="Mass" Name="Mass_myRes">
			<Value>2.0</Value>
			<Error>8.0E-04</Error>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="1"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="+1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Type="Width" Name="Width_myRes">
				<Value>1.0</Value>
        <Min>0.1</Min>
        <Max>1.0</Max>
        <Fix>false</Fix>
			</Parameter>
			<Parameter Type="MesonRadius" Name="Radius_myRes">
				<Value>2.5</Value>
				<Fix>true</Fix>
			</Parameter>
		</DecayInfo>
	</Particle>
</ParticleList>
'''

############################# Start Fit ##############

print("     Create RunManager")

run = RunManager()

  
print("     Create Particle List")

partL = PartList()
ReadParticles(partL, GetDefaultParticles())
ReadParticles(partL, myParticles)


print("     Create Kinematics")

kin = HelicityKinematics(partL, GetInitialState(443), GetFinalState(22, 111, 111))


print("     Generate Phasespace")

gen = RootGenerator(partL, kin)
run.SetGenerator(gen)
phspSample = Data()
run.SetPhspSample(phspSample);
run.GeneratePhsp(100000)


print("     Create Amplitude")

intens = GetIncoherentIntensity(amplitudeModel, partL, kin, phspSample)
run.SetAmplitude(intens)


print("     Generate Data")

sample = Data()
run.SetData(sample)
run.Generate(kin, 1000)


print("     Fit model to data")

fitPar = ParameterList()
intens.GetParameters(fitPar)
setErrorOnParameterList(fitPar, 0.05, False);

esti = MinLogLH(kin, intens, sample, phspSample, phspSample, 0, 0)
esti.UseFunctionTree(True)

minuitif = MinuitIF(esti, fitPar)
minuitif.SetHesse(True)
run.SetOptimizer(minuitif)

result = run.Fit(fitPar)

#fitFracs = CalculateFitFractions(kin, intens, phspPoints, fitComponents);
#CalcFractionError(fitPar, result.GetCovarianceMatrix(), fitFracs, kin,
#                           intens, phspPoints, 100, fitComponents);

result.Print()


print("     Save results")

saveResults("PyDalitzFit-fitResult.xml", result)
saveModel("PyDalitzFit-Model.xml", partL, fitPar, intens)


print("     TODO: Plot ?")

############################# Complete C++ Fit ##############

#print("     test ComPWA c++ fit")

#fit = PythonFit()

#fit.StartFit()

exit()
 