""" sample script for the testing purposes using the decay
    JPsi -> gamma pi0 pi0
"""
import logging

from expertsystem.ui.system_control import StateTransitionManager

from expertsystem.amplitude.helicitydecay import (
    HelicityDecayAmplitudeGeneratorXML)

logging.basicConfig(level=logging.INFO)

# initialize the graph edges (intial and final state)
initial_state = [("D0", [0])]
final_state = [("K_S0", [0]), ("K+", [0]), ("K-", [0])]

tbd_manager = StateTransitionManager(initial_state, final_state,
                                     ['f0', 'a0', 'phi', 'a2(1320)-'])
tbd_manager.number_of_threads = 4

graph_interaction_settings_groups = tbd_manager.prepare_graphs()
(solutions, violated_rules) = tbd_manager.find_solutions(
    graph_interaction_settings_groups)

print("found " + str(len(solutions)) + " solutions!")

# print intermediate state particle names
for g in solutions:
    print(g.edge_props[1]['@Name'])

xml_generator = HelicityDecayAmplitudeGeneratorXML()
xml_generator.generate(solutions)
xml_generator.write_to_file('output2.xml')
