#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pycompwa as pwa

pwa.Logging("out.log", "trace")

# Fill particle list
partL = pwa.PartList()
pwa.read_particles(partL, pwa.default_particles())
partLTrue = pwa.PartList()
pwa.read_particles(partLTrue, pwa.default_particles())
with open('model.xml', 'r') as content_file:
    customPart = content_file.read()
    pwa.read_particles(partL, customPart)
    pwa.read_particles(partLTrue, customPart)

# Create kinematics
kin = pwa.HelicityKinematics(partL, [421], [310, 321, -321])
kin.set_phsp_volume(0.541493)
kinTrue = pwa.HelicityKinematics(partLTrue, [421], [310, 321, -321])
kinTrue.set_phsp_volume(0.541493)

# Generate phase space sample
gen = pwa.RootGenerator(partL, kin, 12345)
phspSample = pwa.Data()
pwa.generate_phsp(100000, gen, phspSample)

# Create Amplitude
with open('model.xml', 'r') as content_file:
    amplitudeModel = content_file.read()
    intens = pwa.incoherent_intensity(amplitudeModel, partL, kin, phspSample,
                                      phspSample)
    intensTrue = pwa.incoherent_intensity(amplitudeModel, partLTrue, kinTrue,
                                          phspSample, phspSample)

kin.print_sub_systems()

# Generate Data
sample = pwa.Data()
# Phase space samples are created on the fly
pwa.generate(1000, kinTrue, gen, intensTrue, sample)
# pwa.generate(1000, kin, gen, intens, sample, phspSample, phspSample)

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

components = [["phi(1020)", "D0toKSK+K-"], ["a0(980)0", "D0toKSK+K-"],
              ["a0(980)+", "D0toKSK+K-"], ["a2(1320)-", "D0toKSK+K-"]]
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
rootpl.add_component(components[0][0], components[0][1], "phi_1020")
rootpl.add_component(components[1][0], components[1][1], "a0_980_0")
rootpl.add_component(components[2][0], components[2][1], "a0_980_plus")
rootpl.add_component(components[3][0], components[3][1], "a2_1320_minus")
rootpl.add_component("BkgD0toKSK+K-", "BkgD0toKSK+K-", "Background")
rootpl.write("", "plot.root", "RECREATE")

exit()
