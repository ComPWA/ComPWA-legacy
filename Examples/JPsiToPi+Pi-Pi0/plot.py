from Plotting.plot import make_dalitz_plot, make_angular_distribution_plot
from Plotting.ROOT.rootplotdatareader import open_compwa_plot_data, load_root


plot_data = open_compwa_plot_data("plot.root")

invariant_masses = plot_data.data.filter(regex='mSq')

#print(invariant_masses[[invariant_masses.columns[0]]])
#print(invariant_masses[[invariant_masses.columns[1]]])
#make_dalitz_plot(invariant_masses[[invariant_masses.columns[0]]], invariant_masses[[invariant_masses.columns[1]]])

make_angular_distribution_plot(plot_data.data[['cosTheta_34_2']], alpha=0.7, bins=40)
make_angular_distribution_plot(plot_data.data[['phi_34_2']], alpha=0.7, bins=40)
make_angular_distribution_plot(plot_data.data[['cosTheta_3_4_vs_2']], alpha=0.7, bins=40)
make_angular_distribution_plot(plot_data.data[['phi_3_4_vs_2']], alpha=0.7, bins=40)
