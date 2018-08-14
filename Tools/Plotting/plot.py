import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import re

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
        counts, bin_edges = np.histogram(
            data_values, kwargs['bins'], weights=data_weights.values)
        bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.
        err = np.sqrt(counts)
        plt.errorbar(bin_centres, counts, yerr=err, fmt='o', label='data')
        axis = plt.gca()
        df = pd.DataFrame({'fit': plot_data.fit_result_data[col_name]})
        df.plot.hist(df, ax=axis, weights=fit_result_weights.values,
                     histtype='step', legend=False, **kwargs)

        xtitle = convert_helicity_column_name_to_title(col_name, plot_data)
        axis.set_xlabel(xtitle)
        axis.set_ylabel('')
        axis.legend()
        plt.tight_layout()
        plt.savefig(replace_particle_ids_with_name(col_name, plot_data)+'.png', bbox_inches='tight')
        #plt.show()


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


def make_dalitz_plot(msq1, msq2, weights=None, labels=[]):
    #invariant_masses = plot_data.data.filter(regex='mSq')
    # print(invariant_masses[[invariant_masses.columns[0]]])
    # print(invariant_masses[[invariant_masses.columns[1]]])
    # make_dalitz_plot(invariant_masses[[invariant_masses.columns[0]]],
    # invariant_masses[[invariant_masses.columns[1]]])

    # plt.hist2d(msq1, msq2, bins=[data_hist.bins_x, data_hist.bins_y],
    #                                      norm=None, weights=weights,
    #                                       range=[data_hist.x_range,data_hist.y_range],
    #                                       alpha=0.8, cmap=plt.get_cmap('viridis'), cmin=0.000001)
    plt.hist2d(msq1, msq2, norm=None, weights=weights, bins=[80, 80],
               alpha=0.8, cmap=plt.get_cmap('viridis'), cmin=0.000001)
    #ax = plt.gca()
    # ax.set_xlabel('$\\theta^{\mathrm{'+label+'}}_x\,/\\mathrm{mrad}$')
    # ax.set_xlabel('$M_1^2$')
    # ax.set_ylabel('$M_2^2$')
    #cbar = plt.colorbar(img, pad=0.0)
    # cbar.set_label(r'$\frac{\mathrm{Model}-\mathrm{Data}}{\mathrm{Error_{Data}}}$')
    plt.colorbar(label='events', pad=0.0)

    plt.xlabel('$M_1^2$')
    plt.ylabel('$M_2^2$')

    # plt.legend()
    # plt.grid(True)

    #numberstring = '{:.3f}'.format(mean(z_err))
    # plt.figtext(0.55, 0.86, '$\mu(L^{\mathrm{err}}_{\mathrm{fit}}) = '+numberstring+'\%$',
    #            verticalalignment='bottom', horizontalalignment='left',
    #            # transform=plt.transAxes,
    #            color='black', fontsize=18)

    plt.show()
    plt.tight_layout()
    plt.savefig("dalitz-plot.png", bbox_inches='tight')
