#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pycompwa as pwa

pwa.Logging("fit.log", "debug")

# Fill particle list
partL = pwa.PartList()
pwa.read_particles(partL, pwa.default_particles())
with open('model.xml', 'r') as content_file:
    customPart = content_file.read()
    pwa.read_particles(partL, customPart)

# Create kinematics
kin = pwa.HelicityKinematics(partL, pwa.initial_state(421),
                             pwa.final_state(310, 321, -321))

# Generate phase space sample
gen = pwa.RootGenerator(partL, kin)
phspSample = pwa.Data()
pwa.generate_phsp(100000, gen, phspSample)

# Create Amplitude
with open('model.xml', 'r') as content_file:
    amplitudeModel = content_file.read()
    intens = pwa.incoherent_intensity(amplitudeModel, partL, kin, phspSample)

# Generate Data
sample = pwa.Data()
pwa.generate(300, kin, gen, intens, sample, phspSample, phspSample)

# Fit model to data
fitPar = pwa.ParameterList()
intens.parameters(fitPar)
pwa.set_parameter_error(fitPar, 0.05, False)

esti = pwa.MinLogLH(kin, intens, sample, phspSample, phspSample, 0, 0)
esti.enable_function_tree(True)
pwa.print_function_tree(esti)

minuitif = pwa.MinuitIF(esti, fitPar)
minuitif.enable_hesse(True)

result = minuitif.minimize(fitPar)

components = [["a0(980)0", "D0toKSK+K-_inc"], ["phi(1020)", "D0toKSK+K-_inc"],
              ["a0(980)+", "D0toKSK+K-"],["a2(1320)-", "D0toKSK+K-"]]
fitFracs = pwa.fit_fractions(kin, intens, phspSample, components)
# pwa.fit_fractions_error(fitPar, result, fitFracs, intens, kin, phspSample, 100)
result.set_fit_fractions(fitFracs)

result.print()

# Save results
pwa.save_results("fitResult.xml", result)
pwa.save_model("fitModel.xml", partL, fitPar, intens)

# Get results
# mydatapoints = pwa.DataPoints(sample, kin)

# Plotting
rootpl = pwa.RootPlot(kin)
rootpl.set_data(sample)
rootpl.set_phsp_sample(phspSample)
rootpl.set_intensity(intens)
rootpl.add_component(components[0][0], "a0_980_0")
rootpl.add_component(components[1][0], "phi_1020")
rootpl.write("", "plot.root", "RECREATE")

exit()
