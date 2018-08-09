import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap

#font = {'family': 'sans-serif',
#        'serif': ['Helvetica'],
#        'weight': 'normal',
#        'size': 18}
#
#rc('text', usetex=True)
#params = {'text.latex.preamble':
#          [r'\usepackage[alsoload=hep,binary-units=true]{siunitx}',
#           r'\usepackage{amsmath}', r'\DeclareSIUnit\c{c}',
#           r'\newunit[per=slash]{\GeVsqcsq}{\giga\eV\tothe{2}\per\c\tothe{2}}'
#           ]}
#plt.rcParams.update(params)
#rc('font', **font)

custom_map = LinearSegmentedColormap.from_list(name='custom_div_cmap',
                                               colors=['b', 'g', 'r'],
                                               N=50)


class PlotData:
    def __init__(self):
        self.data = []
        self.fit_result_data = []
        self.data_integral = 1.0
        self.weighted_phsp_integral = 1.0
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


def make_angular_distribution_plot(dataframe, **kwargs):
    dataframe.plot.hist(**kwargs)
    plt.tight_layout()
    plt.savefig(dataframe.columns.values[0]+'.png', bbox_inches='tight')
    plt.show()


def make_dalitz_plot(msq1, msq2, weights=None, labels=[]):
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
