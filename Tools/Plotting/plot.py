import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from itertools import combinations
import pandas as pd
import re
from math import cos

# font = {'family': 'sans-serif',
#        'serif': ['Helvetica'],
#        'weight': 'normal',
#        'size': 18}
#
#rc('text', usetex=True)
# params = {'text.latex.preamble':
#          [r'\usepackage[alsoload=hep,binary-units=true]{siunitx}',
#           r'\usepackage{amsmath}', r'\DeclareSIUnit\c{c}',
#           r'\newunit[per=slash]{\GeVsqcsq}{\giga\eV\tothe{2}\per\c\tothe{2}}'
#           ]}
# plt.rcParams.update(params)
#rc('font', **font)

custom_map = LinearSegmentedColormap.from_list(name='custom_div_cmap',
                                               colors=['b', 'g', 'r'],
                                               N=50)


class PlotData:
    def __init__(self):
        self.data = None
        self.fit_result_data = None
        self.particle_id_to_name_mapping = {}


def make_difference_plot_1d():
    # if we have a plus minus range of values lets adjust the color range symmetrically
    # if minval < 0 and maxval > 0:
    #    absmax = max(-minval, maxval)
    #    minval = -absmax
    #    maxval = absmax
    pass


def make_difference_plot_2d():
    pass


def make_comparison_plot_1d(plot_data, column_names, **kwargs):
    for col_name in column_names:
        plt.clf()
        data_weights = (plot_data.data.event_weight *
                        plot_data.data.event_efficiency)
        fit_result_weights = (plot_data.fit_result_data.intensity *
                              plot_data.fit_result_data.event_weight *
                              plot_data.fit_result_data.event_efficiency)

        data_values = plot_data.data[col_name].values
        name = col_name
        if 'theta' in col_name:
            name = 'cos' + col_name
            data_values = [cos(x) for x in data_values]
        counts, bin_edges = np.histogram(
            data_values, kwargs['bins'], weights=data_weights.values)
        bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.
        err = np.sqrt(counts)
        plt.errorbar(bin_centres, counts, yerr=err, fmt='o', label='data')
        axis = plt.gca()
        df = pd.DataFrame({'fit': plot_data.fit_result_data[col_name]})
        if 'theta' in col_name:
            df = df.apply(cos, 'columns')
        df.plot.hist(ax=axis, weights=fit_result_weights.values,
                     histtype='step', legend=False, **kwargs)

        xtitle = convert_helicity_column_name_to_title(name, plot_data)
        axis.set_xlabel(xtitle)
        axis.set_ylabel('')
        axis.legend()
        plt.tight_layout()
        plt.savefig(replace_particle_ids_with_name(
            name, plot_data)+'.png', bbox_inches='tight')
        # plt.show()


def replace_particle_ids_with_name(column_name, plot_data):
    replaced_string = column_name
    for index, name in plot_data.particle_id_to_name_mapping.items():
        replaced_string = replaced_string.replace(str(index), name)
    return replaced_string


def convert_helicity_column_name_to_title(column_name, plot_data):
    title = replace_particle_ids_with_name(column_name, plot_data)

    # fix sub and super scripting
    m = re.search('([^_]+)_([^_]+)_([^_]+)(_vs)?_?([^_]+)?_?([^_]+)?', title)
    if m:
        title = m.groups()[0] + \
            '_{' + m.groups()[1] + ',' + m.groups()[2] + '}'
        if m.groups()[4] or m.groups()[5]:
            recoil_string = ''
            if m.groups()[4]:
                recoil_string += m.groups()[4]
            if m.groups()[5]:
                if recoil_string != '':
                    recoil_string += ','
                recoil_string += m.groups()[5]
            title += '^{' + recoil_string + '}'

    # replace mathematical function
    math_functions = ['cos', 'sin', 'tan', 'exp', 'ln', 'log']
    for x in math_functions:
        if x in title:
            title = title.replace(x, x+'(')
            title += ')'

    # TODO: fix the sub and superscripting of the particle charge or spin
    # ex: pi+ -> pi^+ or pi0->pi^0 or f0 -> f_0

    # at latex \ for greek letters
    greek_letters = ['alpha', 'Alpha', 'beta', 'Beta', 'gamma', 'Gamma',
                     'theta', 'Theta', 'phi', 'Phi', 'pi', 'Pi',
                     'psi', 'Psi', 'rho', 'Rho']
    for x in greek_letters:
        if x in title:
            title = title.replace(x, '\\'+x)

    # and finally wrap everything in a math env
    title = '$' + title + '$'

    return title


def make_dalitz_plots(plot_data, **kwargs):
    invariant_mass_names = [x for x in list(plot_data.data.columns.values)
                            if 'mSq' in x and '_vs_' in x]

    data_weights = (plot_data.data.event_weight *
                    plot_data.data.event_efficiency)
    fit_result_weights = (plot_data.fit_result_data.intensity *
                          plot_data.fit_result_data.event_weight *
                          plot_data.fit_result_data.event_efficiency)

    for [im1, im2] in combinations(invariant_mass_names, 2):
        plt.clf()

        msq1 = plot_data.data[im1].values
        msq2 = plot_data.data[im2].values

        plt.hist2d(msq1, msq2, norm=None, weights=data_weights,
                   cmap=plt.get_cmap('viridis'), cmin=0.000001, **kwargs)

        axis = plt.gca()
        xtitle = convert_helicity_column_name_to_title(im1, plot_data)
        axis.set_xlabel(xtitle)
        ytitle = convert_helicity_column_name_to_title(im2, plot_data)
        axis.set_ylabel(ytitle)
        plt.colorbar(label='events', pad=0.0)
        axis.legend()
        plt.tight_layout()
        plt.savefig(replace_particle_ids_with_name(im1, plot_data)
                    + '_vs_' + replace_particle_ids_with_name(im2, plot_data)
                    + '.png', bbox_inches='tight')
        # plt.show()
