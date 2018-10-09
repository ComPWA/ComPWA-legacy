from Plotting.plot import make_dalitz_plots, make_comparison_plot_1d
from Plotting.ROOT.rootplotdatareader import open_compwa_plot_data

plot_data = open_compwa_plot_data("plot.root")

data_variables = list(plot_data.data.columns.values)
print("found data variables:", data_variables)
plot_variables = []
for x in data_variables:
    if 'mSq' in x or 'theta' in x or 'phi' in x:
        plot_variables.append(x)

make_comparison_plot_1d(plot_data, plot_variables, alpha=0.7, bins=40)

make_dalitz_plots(plot_data, alpha=0.7, bins=[50, 50])
