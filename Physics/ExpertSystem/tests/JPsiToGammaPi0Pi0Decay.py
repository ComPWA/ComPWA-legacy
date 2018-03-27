""" sample script for the testing purposes using the decay
    JPsi -> gamma pi0 pi0
"""
import logging

from core.ui.decay_manager import TwoBodyDecayManager

#logging.basicConfig(level=logging.DEBUG)

# initialize the graph edges (initial and final state)
initial_state = [("J/psi", [-1, 1])]
final_state = [("gamma", [-1, 1]), ("pi0", [0]), ("pi0", [0])]

tbd_manager = TwoBodyDecayManager(initial_state, final_state)

graph_node_setting_pairs = tbd_manager.prepare_graphs()
solutions = tbd_manager.find_solutions(graph_node_setting_pairs)

print("found " + str(len(solutions)) + " solutions!")

print("intermediate states:")
for g in solutions:
    print(g.edge_props[1]['@Name'])
