#!/usr/bin/env python3

import pycompwa as pwa
import logging
from math import cos

from expertsystem.ui.system_control import (
    StateTransitionManager, InteractionTypes)

from expertsystem.amplitude.helicitydecay import (
    HelicityDecayAmplitudeGeneratorXML)

from expertsystem.state.particle import get_xml_label, XMLLabelConstants

pwa.Logging("out.log", "trace")

logging.basicConfig(level=logging.INFO)


def generate_model_xml():
    # initialize the graph edges (initial and final state)
    initial_state = [("D1(2420)0", [1])]
    final_state = [("D0", [0]), ("pi-", [0]), ("pi+", [0])]

    tbd_manager = StateTransitionManager(initial_state, final_state,
                                         ['D*'])
    tbd_manager.number_of_threads = 1
    tbd_manager.set_allowed_interaction_types(
        [InteractionTypes.Strong])
    graph_interaction_settings_groups = tbd_manager.prepare_graphs()

    (solutions, violated_rules) = tbd_manager.find_solutions(
        graph_interaction_settings_groups)

    print("found " + str(len(solutions)) + " solutions!")

    print("intermediate states:")
    decinfo_label = get_xml_label(XMLLabelConstants.DecayInfo)
    for g in solutions:
        print(g.edge_props[1]['@Name'])
        for edge_props in g.edge_props.values():
            if decinfo_label in edge_props:
                del edge_props[decinfo_label]
                edge_props[decinfo_label] = {
                    get_xml_label(XMLLabelConstants.Type): "nonResonant"}

    xml_generator = HelicityDecayAmplitudeGeneratorXML()
    xml_generator.generate(solutions)
    xml_generator.write_to_file('model.xml')


def test_angular_distributions(make_plots=False):
    generate_data_samples("model.xml", "plot.root")
    # In this example model the magnitude of A_00 = 0.5 and of A_10=A_-10=1
    # x = cos(theta) distribution from D1 decay should be 1.25 + 0.75*x^2
    # x = cos(theta') distribution from D* decay should be 1 - 0.75*x^2
    # dphi = phi - phi' distribution should be 1 - 1/2.25*cos(2*dphi)
    tuples = [(['theta_34_2'], {'number_of_bins': 200},
               lambda x: 1.25+0.75*x*x),
              (['theta_3_4_vs_2'], {'number_of_bins': 200},
               lambda x: 1-0.75*x*x),
              (['phi_34_2'],
               {'number_of_bins': 200, 'second_column_names': ['phi_3_4_vs_2'],
                'binary_operator': lambda x, y: x-y},
               lambda x: 1-1/2.25*cos(2*x))]
    compare_data_samples_and_theory("plot.root", tuples, make_plots)


def generate_data_samples(model_filename, output_filename):
    # Fill particle list
    partL = pwa.PartList()
    pwa.read_particles(partL, pwa.default_particles())
    partLTrue = pwa.PartList()
    pwa.read_particles(partLTrue, pwa.default_particles())
    with open(model_filename, 'r') as content_file:
        customPart = content_file.read()
        pwa.read_particles(partL, customPart)
        pwa.read_particles(partLTrue, customPart)

    # Create kinematics
    kinTrue = pwa.HelicityKinematics(partLTrue, model_filename)
    # kinTrue.set_phsp_volume(0.541493)

    # Generate phase space sample
    gen = pwa.RootGenerator(partL, kinTrue, 12345)
    phspSample = pwa.generate_phsp(100000, gen)

    # Create Amplitude
    with open(model_filename, 'r') as content_file:
        amplitudeModel = content_file.read()
        intensTrue = pwa.incoherent_intensity(amplitudeModel, partLTrue,
                                              kinTrue, phspSample,
                                              phspSample)

    # Generate Data
    sample = pwa.generate(40000, kinTrue, gen, intensTrue)

    # Plotting
    kinTrue.create_all_subsystems()
    rootpl = pwa.RootPlotData(kinTrue, intensTrue)
    rootpl.set_data(sample)
    rootpl.set_phsp_mc(phspSample)

    rootpl.write("tree", output_filename, "RECREATE")


def compare_data_samples_and_theory(input_rootfile,
                                    distribution_test_tuples,
                                    make_plots):
    from Plotting.plot import (
        make_binned_distributions, chisquare_test, plot_distributions_1d,
        function_to_histogram, scale_to_other_histogram,
        convert_helicity_column_name_to_title
    )
    from Plotting.ROOT.rootplotdatareader import open_compwa_plot_data

    plot_data = open_compwa_plot_data(input_rootfile)

    #data_variables = list(plot_data.data.columns.values)
    #print("found data variables:", data_variables)

    for var_names, kwargs, func in distribution_test_tuples:
        binned_dists = make_binned_distributions(
            plot_data, var_names, **kwargs)
        for var_name, dists in binned_dists.items():
            data_hist = dists['data']
            if make_plots:
                function_hist = function_to_histogram(func, data_hist[0])
                function_hist = scale_to_other_histogram(
                    function_hist, data_hist)

                hist_bundle = {'data': data_hist,
                               'theory': function_hist + ({'fmt': '-'},)
                               }
                xtitle = convert_helicity_column_name_to_title(
                    var_name, plot_data)
                plot_distributions_1d(hist_bundle, var_name, xtitle=xtitle)

            chisquare_value = chisquare_test(data_hist[:-1], func)
            assert(abs(1.0 - chisquare_value) < 0.1)


if __name__ == '__main__':
    test_angular_distributions(make_plots=True)
