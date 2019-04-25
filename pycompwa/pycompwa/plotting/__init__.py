from . import rootplotdatareader
import logging
from math import cos, pi, sqrt
import re
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt


class PlotData:
    def __init__(self, data_record=None,
                 fit_result_data_record=None):
        self.data = data_record
        self.fit_result_data = fit_result_data_record
        self.particle_id_to_name_mapping = {}


def create_nprecord(column_names, data_array):
    if isinstance(data_array, list):
        data_array = np.asarray(data_array)
    if column_names and isinstance(data_array, np.ndarray):
        if len(data_array) == len(column_names):
            return np.rec.fromarrays(data_array, names=column_names)
        else:
            if len(data_array.T) == len(column_names):
                return np.rec.fromarrays(
                    data_array.T, names=column_names)
            else:
                raise ValueError("Data columns and column names mismatch!")


def correct_phi_range(phi):
    while True:
        if phi > pi:
            phi -= 2*pi
        elif phi < -pi:
            phi += 2*pi
        else:
            return phi


def chisquare_test(histogram, func):
    bin_centers, bin_contents, bin_errors = histogram
    function_hist = function_to_histogram(func, bin_centers)
    function_hist = scale_to_other_histogram(function_hist, histogram)
    dof = len(bin_centers)-1
    redchi2 = chisquare(bin_contents, bin_errors,
                        function_hist[1])/dof
    logging.info("chisquare/dof: " + str(redchi2) + " +- " + str(sqrt(2/dof)))

    return redchi2, sqrt(2/dof)


def scale_to_other_histogram(histogram, histogram_reference):
    bin_centers, bin_contents, bin_errors = histogram
    normalization = sum(histogram_reference[1]) / sum(bin_contents)
    logging.info("calculated normalization:" + str(normalization))
    new_bin_contents = [normalization*x for x in bin_contents]
    return (bin_centers, new_bin_contents, np.sqrt(bin_errors))


def function_to_histogram(func, bin_centers):
    import scipy.integrate as integrate
    half_bin_width = (bin_centers[1]-bin_centers[0])/2
    return (
        bin_centers,
        [integrate.quad(
            func, x-half_bin_width, x+half_bin_width)[0] for x in bin_centers],
        []
    )


def chisquare(values, errors, expected):
    return sum([((x[0]-x[1])/x[2])**2
                for x in zip(values, expected, errors)])


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


def make_binned_distributions(plot_data, column_names,
                              number_of_bins=50,
                              use_bin_centers=True,
                              nofit=True,
                              binary_operator=None,
                              second_column_names=None,
                              use_theta=False,
                              **kwargs):
    binned_distributions = {}
    data_weights = (plot_data.data.weight)

    fit_result_weights = np.array([])
    if not nofit and not plot_data.fit_result_data.empty:
        fit_result_weights = (plot_data.fit_result_data.intensity *
                              plot_data.fit_result_data.weight)

    if second_column_names:
        if len(second_column_names) != len(column_names):
            raise ValueError("Number of second column names not ")
        if not binary_operator:
            raise ValueError("binary_operator is not set!")
    else:
        second_column_names = [None for x in range(0, len(column_names))]
    for col_name1, col_name2 in zip(column_names, second_column_names):
        if col_name2:
            data_values = [binary_operator(x, y)
                           for x, y in zip(plot_data.data[col_name1],
                                           plot_data.data[col_name2])]
            if "phi" in col_name1:
                data_values = [correct_phi_range(x) for x in data_values]
            name = col_name1 + "OP" + col_name2
        else:
            data_values = plot_data.data[col_name1]
            name = col_name1
            if not use_theta and 'theta' in col_name1:
                name = 'cos' + col_name1
                data_values = [cos(x) for x in data_values]

        bin_edges, bin_content, errs = make_histogram(
            data_values, data_weights, number_of_bins)
        if use_bin_centers:
            bin_values = (bin_edges[:-1] + bin_edges[1:])/2.
        else:
            bin_values = bin_edges
        binned_distributions[name] = {
            'data': (bin_values, bin_content, errs, {'fmt': 'o'})
        }

        if fit_result_weights.size > 0:
            if col_name2:
                fit_data_values = [
                    binary_operator(x, y)
                    for x, y in zip(
                        plot_data.fit_result_data[col_name1],
                        plot_data.fit_result_data[col_name2])]
                if "phi" in col_name1:
                    fit_data_values = [correct_phi_range(
                        x) for x in fit_data_values]
            else:
                fit_data_values = plot_data.fit_result_data[col_name1]
                if not use_theta and 'theta' in col_name1:
                    fit_data_values = [cos(x) for x in fit_data_values]

            fit_bin_edges, fit_bin_content, fit_errs = make_histogram(
                fit_data_values, fit_result_weights, bin_edges)
            binned_distributions[name]['fit'] = (
                bin_values, fit_bin_content, fit_errs, {})

    return binned_distributions


def make_histogram(values, weights, bins):
    bin_content, bin_edges = np.histogram(
        values, bins=bins, weights=weights)
    errs = np.sqrt(bin_content)

    return (bin_edges, bin_content, errs)


def make_comparison_plot_1d(plot_data, column_names, **kwargs):
    for var_name, distributions in make_binned_distributions(
            plot_data, column_names, **kwargs).items():
        xtitle = convert_helicity_column_name_to_title(var_name, plot_data)
        real_var_name = replace_particle_ids_with_name(var_name, plot_data)
        plot_distributions_1d(distributions, real_var_name, xtitle=xtitle)


def plot_distributions_1d(distributions, var_name, **kwargs):
    plt.clf()
    for name, histogram in distributions.items():
        yerrors = None
        if histogram[2]:
            yerrors = histogram[2]
        plt.errorbar(histogram[0], histogram[1], yerr=yerrors,
                     label=name, **(histogram[3]))

    if plt.ylim()[0] > 0.0:
        plt.ylim(bottom=0.0)
    axis = plt.gca()
    if 'xtitle' in kwargs:
        axis.set_xlabel(kwargs['xtitle'])
    else:
        axis.set_xlabel(var_name)
    axis.set_ylabel('')
    axis.legend()
    plt.tight_layout()
    plt.savefig(var_name+'.png', bbox_inches='tight')


def make_dalitz_plots(plot_data, var_names, **kwargs):
    if var_names:
        invariant_mass_names = var_names
    else:
        invariant_mass_names = [x for x in list(plot_data.data.dtype.names)
                                if 'mSq' in x and '_vs_' in x]

    data_weights = (plot_data.data.weight)

    fit_result_weights = None
    if plot_data.fit_result_data is not None:
        fit_result_weights = (plot_data.fit_result_data.intensity *
                              plot_data.fit_result_data.weight)
        rescale_factor = sum(data_weights)/sum(fit_result_weights)
        fit_result_weights *= rescale_factor

    mass_combinations = list(combinations(invariant_mass_names, 2))

    if fit_result_weights is not None:
        fig, axs = plt.subplots(len(mass_combinations), 2, squeeze=False)
    else:
        fig, axs = plt.subplots(len(mass_combinations), 1, squeeze=False)

    for i, [im1, im2] in enumerate(mass_combinations):
        msq1 = plot_data.data[im1]
        msq2 = plot_data.data[im2]

        if 'cmin' not in kwargs:
            kwargs['cmin'] = 0.00001

        img = axs[i, 0].hist2d(msq1, msq2, norm=None,
                               weights=data_weights,
                               cmap=plt.get_cmap('viridis'),
                               **kwargs)

        xtitle = convert_helicity_column_name_to_title(im1, plot_data)
        axs[i, 0].set_xlabel(xtitle)
        ytitle = convert_helicity_column_name_to_title(im2, plot_data)
        axs[i, 0].set_ylabel(ytitle)
        plt.colorbar(img[3], ax=axs[i, 0],
                     label='events', orientation='vertical', pad=0.0)
        # axis.legend()
        # axs[0, i].savefig(replace_particle_ids_with_name(im1, plot_data)
        #            + '_vs_' + replace_particle_ids_with_name(im2, plot_data)
        #            + '.png', bbox_inches='tight')

        if fit_result_weights is not None:
            msq1 = plot_data.fit_result_data[im1]
            msq2 = plot_data.fit_result_data[im2]
            img = axs[i, 1].hist2d(msq1, msq2, norm=None,
                                   weights=fit_result_weights,
                                   cmap=plt.get_cmap('viridis'),
                                   **kwargs)

            xtitle = convert_helicity_column_name_to_title(im1, plot_data)
            axs[i, 1].set_xlabel(xtitle)
            ytitle = convert_helicity_column_name_to_title(im2, plot_data)
            axs[i, 1].set_ylabel(ytitle)
            plt.colorbar(
                img[3], ax=axs[i, 1], label='events',
                orientation='vertical', pad=0.0)

    plt.tight_layout()


def make_difference_plot_1d():
    # if we have a plus minus range of values lets adjust the color range
    # symmetrically
    # if minval < 0 and maxval > 0:
    #    absmax = max(-minval, maxval)
    #    minval = -absmax
    #    maxval = absmax
    pass


def make_difference_plot_2d(plot_data, var_names, **kwargs):
    if not isinstance(var_names, (list, tuple)) or not len(var_names) == 2:
        raise ValueError(
            "Incorrent number of variable names! Expecting two variables.")

    data_weights = (plot_data.data.weight)

    if plot_data.fit_result_data is None:
        raise ValueError("Fit result data has to be present!")

    fit_result_weights = (plot_data.fit_result_data.intensity *
                          plot_data.fit_result_data.weight)

    rescale_factor = sum(data_weights)/sum(fit_result_weights)
    fit_result_weights *= rescale_factor

    vars1 = plot_data.data[var_names[0]]
    vars2 = plot_data.data[var_names[1]]

    h, xedges, yedges = np.histogram2d(vars1, vars2, normed=None,
                                       weights=data_weights,
                                       **kwargs)

    hfit, xedges, yedges = np.histogram2d(
        plot_data.fit_result_data[var_names[0]],
        plot_data.fit_result_data[var_names[1]],
        bins=[xedges, yedges], normed=None,
        weights=fit_result_weights)

    plt.clf()
    hdiff = (h-hfit).T

    minweight = hdiff.min()
    maxweight = hdiff.max()
    maxweight = max([abs(minweight), maxweight])

    plt.imshow(hdiff, interpolation='nearest', origin='low',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               cmap=plt.get_cmap('bwr'), vmin=-maxweight, vmax=maxweight)

    plt.colorbar(ax=plt.gca(), label='data-model',
                 orientation='vertical', pad=0.0)
    plt.tight_layout()
