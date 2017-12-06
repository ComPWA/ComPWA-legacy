#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ROOT
from ROOT import gROOT, TCanvas, TF1, TVector3, TTree

from rootpy.interactive import wait

import numpy as np
import matplotlib.pyplot as plt

import PyComPWA as pwa

import time

from histogrammar.tutorial import cmsdata
events = cmsdata.EventIterator()

from histogrammar import *

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
  
print("     Create Particle List")

partL = pwa.PartList()
pwa.ReadParticles(partL, pwa.GetDefaultParticles())
pwa.ReadParticles(partL, myParticles)


print("     Create Kinematics")

kin = pwa.HelicityKinematics(partL, pwa.GetInitialState(443), pwa.GetFinalState(22, 111, 111))


print("     Generate Phasespace")

gen = pwa.RootGenerator(partL, kin)
phspSample = pwa.Data()
pwa.GeneratePhsp(100000, gen, phspSample)


print("     Create Amplitude")

intens = pwa.GetIncoherentIntensity(amplitudeModel, partL, kin, phspSample)


print("     Generate Data")

sample = pwa.Data()
pwa.Generate(1000, kin, gen, intens, sample, phspSample, phspSample)


print("     Fit model to data")

fitPar = pwa.ParameterList()
intens.GetParameters(fitPar)
pwa.setErrorOnParameterList(fitPar, 0.05, False)
print(fitPar.ToString())

esti = pwa.MinLogLH(kin, intens, sample, phspSample, phspSample, 0, 0)
esti.UseFunctionTree(True)

minuitif = pwa.MinuitIF(esti, fitPar)
minuitif.SetHesse(True)

result = minuitif.exec(fitPar)

fitFracs = pwa.CalculateFitFractions(kin, intens, phspSample)
pwa.CalcFractionError(fitPar, result, fitFracs, intens, kin,
                            phspSample, 100)
result.SetFitFractions(fitFracs)  

result.Print()


print("     Save results")

pwa.saveResults("PyDalitzFit-fitResult.xml", result)
pwa.saveModel("PyDalitzFit-Model.xml", partL, fitPar, intens)

print("     Get results")

mydatapoints = pwa.DataPoints(sample, kin)
mynumpydata = np.array(mydatapoints, copy = False)

print("     Plot")

#resultPar = result.GetFinalParameters()

parPlotX = np.arange(0, mynumpydata.shape[0], 1)
parPlotY = mynumpydata[:,0]
plt.plot(parPlotX, parPlotY)

plt.xlabel('Event ID')
plt.ylabel('Event weight')
plt.title('Event weights')
plt.grid(True)
plt.savefig("testPyPlot.png")
plt.show()

#rootpl = pwa.RootPlot(kin)
#rootpl.SetData(sample)
#rootpl.SetPhspSample(phspSample)
#rootpl.SetFitAmp(intens)
#rootpl.Write("PyTestFit", "PyRootPlot.root", "RECREATE")

#resultFile = ROOT.TFile("DalitzFit.root")
#invMassPlot = resultFile.Get("invmass")
#invMassPlot.Draw()
#wait()
#time.sleep(60)

#hist2d = Bin(10, -100, 100, lambda event: event.met.px,
#             value = Bin(10, -100, 100, lambda event: event.met.py))

#for i, event in enumerate(events):
#    if i == 1000: break
#    hist2d.fill(event)

#roothist = hist2d.plot.root("name2", "title")
#roothist.Draw("colz")

exit()
 