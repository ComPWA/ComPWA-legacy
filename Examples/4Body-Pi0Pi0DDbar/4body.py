#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pycompwa as pwa

pwa.Logging("out.log", "trace")

# Fill particle list
partL = pwa.PartList()
pwa.read_particles(partL, pwa.default_particles())
partLTrue = pwa.PartList()
pwa.read_particles(partLTrue, pwa.default_particles())
with open('4body-model.xml', 'r') as content_file:
    customPart = content_file.read()
    pwa.read_particles(partL, customPart)
    pwa.read_particles(partLTrue, customPart)

# Create kinematics
kin = pwa.HelicityKinematics(partL, [22], [111, 111, 421, -421],
                             [0, 0, 0, 4.2])
kin.set_phsp_volume(0.541493)
kinTrue = pwa.HelicityKinematics(partLTrue, [22], [111, 111, 421, -421],
                                 [0, 0, 0, 4.2])
kinTrue.set_phsp_volume(0.541493)

# Generate phase space sample
gen = pwa.RootGenerator(partL, kin, 12345)
phspSample = pwa.Data()
pwa.generate_phsp(100000, gen, phspSample)

# Create Amplitude
with open('4body-model.xml', 'r') as content_file:
    amplitudeModel = content_file.read()
    intens = pwa.incoherent_intensity(amplitudeModel, partL, kin, phspSample,
                                      phspSample)
    intensTrue = pwa.incoherent_intensity(amplitudeModel, partLTrue, kinTrue,
                                          phspSample, phspSample)

kin.print_sub_systems()

# Generate Data
sample = pwa.Data()
# Phase space samples are created on the fly
pwa.generate(50, kinTrue, gen, intensTrue, sample)

# Fit model to data
truePar = pwa.ParameterList()
intensTrue.parameters(truePar)
fitPar = pwa.ParameterList()
intens.parameters(fitPar)
fitPar.set_parameter_error(0.06, False)
pwa.log(fitPar)

esti = pwa.MinLogLH(kin, intens, sample, phspSample, phspSample, 0, 0)
esti.enable_function_tree(True)
esti.log_function_tree()

minuitif = pwa.MinuitIF(esti, fitPar)
minuitif.enable_hesse(True)

result = minuitif.minimize(fitPar)

components = [["Y(1)ToD*0(-1)D*0bar(-1)ToPi0Pi0D0D0bar",
               "Y(1)ToD*0D*0barToPi0Pi0D0bar"]]
fitFracs = pwa.fit_fractions(kin, intens, phspSample, components)
pwa.log(fitFracs)
pwa.fit_fractions_error(fitPar, result, fitFracs, intens, components, kin,
                        phspSample, 20)
result.set_fit_fractions(fitFracs)

result.print()

# Save results
result.write("out-fitResult.xml")
intens.write("out-fitModel.xml")
print(fitPar)
partL.update(fitPar)
partL.write("out-fitParticles.xml")

# Get results
# mydatapoints = pwa.DataPoints(sample, kin)

# Plotting
rootpl = pwa.RootPlot(kin)
rootpl.set_data(sample)
rootpl.set_phsp_sample(phspSample)
rootpl.set_intensity(intens)
rootpl.add_component(components[0][0], components[0][1], "-1_-1")
rootpl.write("", "plot.root", "RECREATE")

exit()
