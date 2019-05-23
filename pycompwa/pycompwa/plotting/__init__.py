import logging
from math import cos, pi, sqrt
import re
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import copy


class PlotData:
    def __init__(self, data_record=np.array([]),
                 fit_result_data_record=np.array([])):
        self.data = data_record
        self.fit_result_data = fit_result_data_record
        self.particle_id_to_name_mapping = {}


class Histogram:
    def __init__(self, dimensions, bin_edges, bin_contents, bin_errors=None,
                 **mpl_kwargs):
        self.dimensions = dimensions
        self.bin_edges = bin_edges
        self.bin_contents = bin_contents
        self.bin_errors = bin_errors
        self.mpl_kwargs = mpl_kwargs


class Dimension:
    def __init__(self, field_name, binary_operator=None,
                 secondary_field_name=None):
        self.field_name = field_name
        self.binary_operator = binary_operator
        self.secondary_field_name = secondary_field_name


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
    from functools import reduce
    function_hist = function_to_histogram(func, histogram)
    function_hist = scale_to_other_histogram(function_hist, histogram)
    dof = reduce((lambda x, y: x * y),
                 np.asarray(histogram.bin_contents).shape)-1

    redchi2 = chisquare(histogram.bin_contents, histogram.bin_errors,
                        function_hist.bin_contents)/dof
    logging.info("chisquare/dof: " + str(redchi2) +
                 " +- " + str(2.0*sqrt(2/dof)))

    return 1.0, redchi2, 2.0*sqrt(2/dof)


def function_to_histogram(func, histogram):
    bin_edges = histogram.bin_edges
    int_ranges = [[x] for x in zip(bin_edges[0][:-1], bin_edges[0][1:])]
    for dim_edges in bin_edges[1:]:
        new_int_ranges = []
        for x in int_ranges:
            temp_int_ranges_row = []
            for new_range in zip(dim_edges[:-1], dim_edges[1:]):
                temprange = copy.copy(x)
                temprange.append(new_range)
                temp_int_ranges_row.append(temprange)
            new_int_ranges.append(temp_int_ranges_row)
        int_ranges = new_int_ranges

    return Histogram(histogram.dimensions, bin_edges,
                     integrate_row_of_bins(func, int_ranges))


def integrate_row_of_bins(func, integration_ranges):
    import scipy.integrate as integrate
    if isinstance(integration_ranges[0][0], tuple):
        return [integrate.nquad(func, x)[0] for x in integration_ranges]
    return [integrate_row_of_bins(func, x) for x in integration_ranges]


def scale_to_other_histogram(histogram, histogram_reference):
    normalization = np.sum(histogram_reference.bin_contents) / \
        np.sum(histogram.bin_contents)
    logging.info("calculated normalization:" + str(normalization))
    new_bin_contents = np.multiply(normalization, histogram.bin_contents)
    new_bin_errors = None
    if histogram.bin_errors:
        new_bin_errors = [np.sqrt(normalization) *
                          x for x in histogram.bin_errors]
    return Histogram(histogram.dimensions, histogram.bin_edges,
                     new_bin_contents, bin_errors=new_bin_errors)


def chisquare(values, errors, expected):
    return np.sum([((x[0]-x[1])/x[2])**2
                   for x in zip(values, expected, errors)])


def replace_particle_ids_with_name(column_name, plot_data):
    replaced_string = column_name
    if plot_data:
        for index, name in plot_data.particle_id_to_name_mapping.items():
            replaced_string = replaced_string.replace(str(index), name)
    return replaced_string


def create_axis_title(dimension, plot_data):
    dimension_label = dimension.field_name
    if dimension.binary_operator and dimension.secondary_field_name:
        dimension_label += '*OP*' + dimension.secondary_field_name
    title = replace_particle_ids_with_name(dimension_label, plot_data)

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


def make_binned_distributions(plot_data, dimensions_list,
                              number_of_bins=50,
                              nofit=False,
                              use_cosine_theta=True,
                              **kwargs):
    binned_distributions = []

    if isinstance(dimensions_list, str):
        dimensions_list = [[Dimension(dimensions_list)]]
    elif isinstance(dimensions_list, Dimension):
        dimensions_list = [[dimensions_list]]
    elif isinstance(dimensions_list, list):
        new_dimensions_list = []
        for x in dimensions_list:
            temp_dim_list = []
            if isinstance(x, list):
                for y in x:
                    temp_dim_list.append([Dimension(x) if isinstance(
                        x, str) else x for x in dimensions_list])
        if not new_dimensions_list:
            new_dimensions_list = [[Dimension(x) if isinstance(
                        x, str) else x for x in dimensions_list]]
        dimensions_list = new_dimensions_list

    data_weights = (plot_data.data.weight)

    fit_result_weights = np.array([])
    if plot_data.fit_result_data.size > 0:
        fit_result_weights = (plot_data.fit_result_data.intensity *
                              plot_data.fit_result_data.weight)
        scaling_factor = sum(data_weights) / sum(fit_result_weights)
        logging.info("scaling fit to data using factor: " +
                     str(scaling_factor))
        fit_result_weights = (fit_result_weights * scaling_factor)

    for dimensions in dimensions_list:
        new_dimensions, data_array = create_data_values(
            plot_data.data, dimensions, use_cosine_theta)
        temp_hist = make_histogram(new_dimensions, data_array, data_weights,
                                   number_of_bins, fmt='o')
        new_distributions = {
            'data': temp_hist
        }

        if fit_result_weights.size > 0:
            new_dimensions, fit_data_array = create_data_values(
                plot_data.fit_result_data, dimensions, use_cosine_theta)
            new_distributions['fit'] = make_histogram(
                new_dimensions, fit_data_array, fit_result_weights,
                bins=temp_hist.bin_edges)

        binned_distributions.append(new_distributions)

    return binned_distributions


def create_data_values(datarecord, dimensions, use_cosine_theta):
    data_array = []
    new_dimensions = []
    for dimension in dimensions:
        col_name = dimension.field_name
        sec_col_name = dimension.secondary_field_name
        if dimension.binary_operator and sec_col_name:
            values1, name1 = change_theta_to_cosine(
                datarecord, col_name, use_cosine_theta)
            values2, name2 = change_theta_to_cosine(
                datarecord, sec_col_name, use_cosine_theta)
            new_values = [dimension.binary_operator(x, y)
                          for x, y in zip(values1, values2)]
            dimension.field_name = name1
            dimension.secondary_field_name = name2
            new_dimensions.append(dimension)
            if "phi" in name1:
                new_values = [correct_phi_range(x) for x in new_values]
            data_array.append(new_values)
        else:
            values1, name1 = change_theta_to_cosine(
                datarecord, col_name, use_cosine_theta)
            dimension.field_name = name1
            new_dimensions.append(dimension)
            data_array.append(values1)

    return (new_dimensions, data_array)


def change_theta_to_cosine(datarecord, column_name, use_cosine_theta):
    if use_cosine_theta:
        if 'theta' in column_name and 'cos' not in column_name:
            return ([cos(x) for x in datarecord[column_name]],
                    'cos'+column_name)
    return (datarecord[column_name], column_name)


def make_histogram(dimensions, values, weights, bins=50, **kwargs):
    bin_content, bin_edges = np.histogramdd(
        values, bins=bins, weights=weights)
    if len(bin_content.shape) == 1:
        errs = [np.sqrt(x) if x > 0 else 1 for x in bin_content]
    elif len(bin_content.shape) == 2:
        errs = [[np.sqrt(x) if x > 0 else 1 for x in row]
                for row in bin_content]
    return Histogram(dimensions, bin_edges, bin_content, errs, **kwargs)


def make_comparison_plot_1d(plot_data, column_names, **kwargs):
    for histograms in make_binned_distributions(
            plot_data, column_names, **kwargs):
        dimension = histograms['data'].dimensions[0]
        xtitle = create_axis_title(dimension, plot_data)
        plot_distributions_1d(histograms, xtitle=xtitle)


def plot_distributions_1d(histograms, use_bin_centers=True,
                          **kwargs):
    plt.clf()
    var_name = ''
    for name, histogram in histograms.items():
        bincenters = histogram.bin_edges
        if use_bin_centers:
            bincenters = 0.5 * \
                (histogram.bin_edges[0][1:]+histogram.bin_edges[0][:-1])
        plt.errorbar(bincenters, histogram.bin_contents,
                     yerr=histogram.bin_errors,
                     label=name, **(histogram.mpl_kwargs))
        if var_name == '':
            var_name = histogram.dimensions[0].field_name

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

    fit_result_weights = np.array([])
    if plot_data.fit_result_data.size > 0:
        fit_result_weights = (plot_data.fit_result_data.intensity *
                              plot_data.fit_result_data.weight)
        rescale_factor = sum(data_weights)/sum(fit_result_weights)
        fit_result_weights *= rescale_factor

    mass_combinations = list(combinations(invariant_mass_names, 2))

    if fit_result_weights.size > 0:
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

        xtitle = create_axis_title(Dimension(im1), plot_data)
        axs[i, 0].set_xlabel(xtitle)
        ytitle = create_axis_title(Dimension(im2), plot_data)
        axs[i, 0].set_ylabel(ytitle)
        plt.colorbar(img[3], ax=axs[i, 0],
                     label='events', orientation='vertical', pad=0.0)

        if fit_result_weights.size > 0:
            msq1 = plot_data.fit_result_data[im1]
            msq2 = plot_data.fit_result_data[im2]
            img = axs[i, 1].hist2d(msq1, msq2, norm=None,
                                   weights=fit_result_weights,
                                   cmap=plt.get_cmap('viridis'),
                                   **kwargs)

            axs[i, 1].set_xlabel(xtitle)
            axs[i, 1].set_ylabel(ytitle)
            plt.colorbar(
                img[3], ax=axs[i, 1], label='events',
                orientation='vertical', pad=0.0)

    plt.tight_layout()


def make_difference_plot_2d(plot_data, var_names, **kwargs):
    if not isinstance(var_names, (list, tuple)) or not len(var_names) == 2:
        raise ValueError(
            "Incorrent number of variable names! Expecting two variables.")

    if plot_data.fit_result_data is None:
        raise ValueError("Fit result data has to be present!")

    dists = make_binned_distributions(plot_data, var_names)

    for dist_pair in dists:
        plot_histogram_difference_2d(dist_pair, **kwargs)


def plot_histogram_difference_2d(histograms, **kwargs):
    keys = list(histograms.keys())
    hdiff = copy.deepcopy(histograms[keys[0]])
    label = keys[0] + '-' + keys[1]
    hdiff.bin_contents = (histograms[keys[0]].bin_contents -
                          histograms[keys[1]].bin_contents)

    plot_distribution_2d(hdiff, True, label, **kwargs)


def plot_distribution_2d(histogram, is_difference=False, zaxis_label=None,
                         **kwargs):
    plt.clf()

    plot_name = histogram.dimensions[0].field_name + \
        '_vs_' + histogram.dimensions[1].field_name

    minweight = histogram.bin_contents.min()
    maxweight = histogram.bin_contents.max()

    colormap = plt.get_cmap('viridis')
    if is_difference:
        plot_name += "_diff"
        maxweight = max([abs(minweight), maxweight])
        minweight = -maxweight
        colormap = plt.get_cmap('bwr')

    xedges = histogram.bin_edges[0]
    yedges = histogram.bin_edges[1]

    plt.imshow(histogram.bin_contents.T, interpolation='nearest', origin='low',
               extent=[xedges[0], xedges[-1], yedges[0],
                       yedges[-1]], aspect='auto',
               cmap=colormap, vmin=minweight, vmax=maxweight)

    axis = plt.gca()
    if 'xtitle' in kwargs:
        axis.set_xlabel(kwargs['xtitle'])
    else:
        axis.set_xlabel(create_axis_title(histogram.dimensions[0], None))
    if 'ytitle' in kwargs:
        axis.set_ylabel(kwargs['ytitle'])
    else:
        axis.set_ylabel(create_axis_title(histogram.dimensions[1], None))

    if not zaxis_label:
        zaxis_label = 'entries'
    plt.colorbar(ax=axis, label=zaxis_label,
                 orientation='vertical', pad=0.0)
    plt.tight_layout()
    plt.savefig(plot_name+'.png', bbox_inches='tight')
