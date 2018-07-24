#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pycompwa as pwa
import logging

from expertsystem.ui.system_control import (
    StateTransitionManager, InteractionTypes)

from expertsystem.amplitude.helicitydecay import (
    HelicityDecayAmplitudeGeneratorXML)

pwa.Logging("out.log", "trace")


logging.basicConfig(level=logging.INFO)

# initialize the graph edges (initial and final state)
initial_state = [("J/psi", [-1, 1])]
final_state = [("pi+", [0]), ("pi-", [0]), ("pi0", [0])]

tbd_manager = StateTransitionManager(initial_state, final_state,
                                     ['rho', 'f'])
#tbd_manager.number_of_threads = 1
tbd_manager.set_interaction_settings(
    [InteractionTypes.EM])
graph_interaction_settings_groups = tbd_manager.prepare_graphs()

(solutions, violated_rules) = tbd_manager.find_solutions(
    graph_interaction_settings_groups)

print("found " + str(len(solutions)) + " solutions!")

print("intermediate states:")
for g in solutions:
    print(g.edge_props[1]['@Name'])

xml_generator = HelicityDecayAmplitudeGeneratorXML()
xml_generator.generate(solutions)
xml_generator.write_to_file('model.xml')



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
kin = pwa.HelicityKinematics(partL, 'model.xml')
#kin.set_phsp_volume(0.541493)
kinTrue = pwa.HelicityKinematics(partLTrue, 'model.xml')
#kinTrue.set_phsp_volume(0.541493)

# Generate phase space sample
gen = pwa.RootGenerator(partL, kin, 12345)
phspSample = pwa.Data()
pwa.generate_phsp(10000, gen, phspSample)

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

#components = [["Y(1)ToD*0(-1)D*0bar(-1)ToPi0Pi0D0D0bar",
#               "Y(1)ToD*0D*0barToPi0Pi0D0bar"]]
#fitFracs = pwa.fit_fractions(kin, intens, phspSample, components)
#pwa.log(fitFracs)
#pwa.fit_fractions_error(fitPar, result, fitFracs, intens, components, kin,
#                        phspSample, 20)
#result.set_fit_fractions(fitFracs)

#result.print()

# Save results
result.write("out-fitResult.xml")
intens.write("out-fitModel.xml")
print(fitPar)
partL.update(fitPar)
partL.write("out-fitParticles.xml")

# Get results
# mydatapoints = pwa.DataPoints(sample, kin)

# Plotting
rootpl = pwa.RootPlotData(kin, intens)
rootpl.set_data(sample)
rootpl.set_phsp_mc(phspSample)
#rootpl.add_component(components[0][0], components[0][1], "-1_-1")
rootpl.write("tree", "plot.root", "RECREATE")

exit()
