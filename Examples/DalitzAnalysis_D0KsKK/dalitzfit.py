import pycompwa.ui as pwa

pwa.Logging("out.log", "trace")

# Create intensity and kinematics
intens, kin = pwa.create_intensity_and_kinematics('model.xml')

# Generate phase space sample
gen = pwa.RootGenerator(
    kin.get_particle_state_transition_kinematics_info(), 12345)
phspSample = pwa.generate_phsp(100000, gen)

kin.print_sub_systems()

# Generate Data
sample = pwa.generate(1000, kin, gen, intens)


# Fit model to data
truePar = pwa.ParameterList()
intens.add_unique_parameters_to(truePar)

fitPar = pwa.ParameterList()
intens.add_unique_parameters_to(fitPar)
fitPar.set_parameter_error(0.06, False)
pwa.log(fitPar)

# convert data to parameter list structure used by the function tree estimator
phspSample.convert_events_to_parameterlist(kin)
sample.convert_events_to_parameterlist(kin)

esti = pwa.create_unbinned_log_likelihood_function_tree_estimator(
    intens, sample, phspSample)

minuitif = pwa.MinuitIF(esti, fitPar)
minuitif.enable_hesse(True)

result = minuitif.minimize(fitPar)

print(fitPar)

# TODO: fix the fit fractions code block as soon as the new implementation is
# available
# components = [["phi(1020)", "D0toKSK+K-"], ["a0(980)0", "D0toKSK+K-"],
#              ["a0(980)+", "D0toKSK+K-"], ["a2(1320)-", "D0toKSK+K-"]]
# fitFracs = pwa.fit_fractions(intens, phspSample, components)
# pwa.log(fitFracs)
# pwa.fit_fractions_error(fitPar, result, fitFracs, intens, components, kin,
#                         phspSample, 20)
# result.set_fit_fractions(fitFracs)

result.log()

# Save results
# TODO: use new implementation, once available

# Plotting
# TODO: fix component part with fit fractions above
pwa.create_rootplotdata("plot.root", kin, sample,
                        phspSample, intens, tfile_option="RECREATE")
# rootpl.add_component(components[0][0], components[0][1], "phi_1020")
# rootpl.add_component(components[1][0], components[1][1], "a0_980_0")
# rootpl.add_component(components[2][0], components[2][1], "a0_980_plus")
# rootpl.add_component(components[3][0], components[3][1], "a2_1320_minus")
# rootpl.add_component("BkgD0toKSK+K-", "BkgD0toKSK+K-", "Background")

exit()
