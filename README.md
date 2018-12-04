[![Build Status](https://travis-ci.org/ComPWA/ComPWA.svg?branch=master)](https://travis-ci.org/ComPWA/ComPWA)
[![Documentation Status](https://readthedocs.org/projects/compwa/badge/?version=latest)](https://compwa.readthedocs.io/en/latest/?badge=latest)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/13697/badge.svg)](https://scan.coverity.com/projects/compwa-compwa)

[![ComPWA Logo](https://github.com/ComPWA/ComPWA/blob/master/doc/images/logo.png)](#)

## About
ComPWA is a project to provide a flexible and modular Partial Wave Analysis framework for various use-cases. It was initiated and is developed by the Panda Collaboration (antiProton ANnihilation at DArmstadt) at Fair (Facility for Antiproton and Ion Research). But ComPWA will not just be used for the Panda physics program, but also for various other experiments to provide a commonly used tool which is stable, efficient and provides comparable results. At the moment there are many PWA-tools on the market, most used just for specific experiments and specific physics cases, some experiments even have multiple tools. But why write the same software again and again? E.g. the model describing physical processes should stay the same independent where and how there was a measurement of the process. Using the actual same implementation of the model does not only save a lot of time, it also ensures that two experiments are comparing the same thing. The same argument holds for optimization-routines and estimation-functions. It might even allow combined fitting of different experiments instead of taking the average of the results!

The natural modularization, following the considerations above, would be to separate into experiment specific information, physics (models, formalisms), estimation how good the model fits the data and optimization of free parameters. The first considerations on this where discussed with experts from different experiments and different technologies where discussed and tested. The result of this process is the first requirement document of the new PWA-Framework.
This sketch illustrates the modular concept: 
[![ComPWA Modules](https://github.com/ComPWA/ComPWA/blob/master/doc/images/compwa_modules.png)](#)

## Available Features :sparkles:
#### Physic Models:
- Helicity Formalism
- Canonical Formalism
#### Data Formats: 
- ROOT
#### Estimators:
- Chi Square
- Unbinned-LogLikelihood
#### Optimizers:
- Miniuit2
- Geneva
#### User Interfaces/Steering:
- C++
- Python

## Example

More detailed examples can be found in the [Examples](https://github.com/ComPWA/ComPWA/tree/master/Examples) directory (especially recommendable are the [jupyter](https://github.com/ComPWA/ComPWA/tree/master/Examples/jupyter) notebooks)!

Before you start make sure your python environment is [correctly](https://github.com/ComPWA/ComPWA/wiki/Installation#setting-up-a-python-virtual-environment) set up.

### Example: J/Psi -> gamma pi0 pi0 decay
In this quickstart example it is shown how to use ComPWA via the python interface. The workflow is:

- (1-4). create a model for the decay first.
- (5.) Afterwards a Monte Carlo data sample (hit & miss) is generated using this model. 
- (6.) Then it is shown how to visualize the data. 
- (7.) Finally a fit on the data sample using the Minuit2 interface is performed. 

Let's go!

`import` the necessary expert system module parts


```python
from expertsystem.ui.system_control import (
    StateTransitionManager, InteractionTypes)
from expertsystem.amplitude.helicitydecay import (
    HelicityDecayAmplitudeGeneratorXML)
from expertsystem.topology.graph import (
    get_intermediate_state_edges)

# just a little function to print the intermediate states
def print_intermediate_states(solutions):
    print("intermediate states:")
    intermediate_states = set()
    for g in solutions:
        edge_id = get_intermediate_state_edges(g)[0]
        intermediate_states.add(g.edge_props[edge_id]['@Name'])
    print(intermediate_states)
```

#### 1. Define problem set

First we define the boundary conditions of our physics problem, such as
- initial state
- final state
- formalism type
- ...

Pass all of that information to the `StateTransitionManager`, which is the main python user interface class of the ComPWA expert system


```python
initial_state = [("J/psi", [-1, 1])]
final_state = [("gamma", [-1, 1]), ("pi0", [0]), ("pi0", [0])]

tbd_manager = StateTransitionManager(initial_state, final_state,
                                     formalism_type='helicity',
                                     topology_building='isobar')
```

#### 2+3. Preparation
Create all topology graphs using the isobar model (two-body decays).

Also initialize the graphs with the initial and final state. Remember that each interaction node defines their own set of conservation laws. The `StateTransitionManager` (STM) defines three interaction types:

| Interaction | Strength |
| --- | --- |
| strong | 60 |
| electromagnetic (EM) | 1 |
| weak | 10^-4 |

Be default all three are used in the preparation stage. `prepare_graphs()` of the STM generates graphs with all possible combinations of interaction nodes. An overall interaction strength is assigned to each graph, and they are grouped according to this strength.


```python
graph_interaction_settings_groups = tbd_manager.prepare_graphs()
```

#### 4. Finding solutions
If you are happy with the automatic settings generated by the StateTransitionManager, just go directly to the solving!


```python
(solutions, violated_rules) = tbd_manager.find_solutions(
        graph_interaction_settings_groups)

print("found " + str(len(solutions)) + " solutions!")
print_intermediate_states(solutions)
```

Ok, now we have a lot of solutions that are actually heavily supressed (involve two weak decays). In general you can modify the dictionary return by `prepare_graphs()` directly.

The STM also comes with a functionality to globally choose the allowed interaction types (`set_allowed_interaction_types()`). Go ahead and **disable** the **EM** and **weak** interaction!


```python
tbd_manager.set_allowed_interaction_types(
    [InteractionTypes.Strong])
graph_interaction_settings_groups = tbd_manager.prepare_graphs()
(solutions, violated_rules) = tbd_manager.find_solutions(
        graph_interaction_settings_groups)
print("found " + str(len(solutions)) + " solutions!")
```

Huh, what happened here? Actually, since a **gamma particle** appears, the expert system knows that there must be **at least one EM interaction** involved. As a consequence no graphs are prepared for this setting!


```python
print(graph_interaction_settings_groups)
```

So let's include the EM interaction...


```python
tbd_manager.set_allowed_interaction_types(
    [InteractionTypes.Strong, InteractionTypes.EM])
graph_interaction_settings_groups = tbd_manager.prepare_graphs()
(solutions, violated_rules) = tbd_manager.find_solutions(
        graph_interaction_settings_groups)

print("found " + str(len(solutions)) + " solutions!")
print_intermediate_states(solutions)
```

Great! Now we selected the solutions that are contribution the strongest. However, be aware that there are more effects that can suppress certain decays. For example branching ratios. In this example **J/Psi** can decay into **pi0 + rho0** or **pi0 + omega**.

| decay | branching ratio |
| --- | --- |
| omega -> gamma+pi0 | 0.0828 |
| rho0 -> gamma+pi0 | 0.0006 |

Unfortunately the **rho0** decays mainly into **pi+pi**, not gamma+pi0. Hence it is suppressed. This information is currently not known to the expert system.
But you can also tell the expert system, which particles you want to allow as intermediate states.


```python
# particles are found by name comparison; so i.e. f2 will find all f2's and f all f's
tbd_manager.allowed_intermediate_particles = ['f']
#tbd_manager.allowed_intermediate_particles = ['f2, f0']

(solutions, violated_rules) = tbd_manager.find_solutions(
        graph_interaction_settings_groups)

print("found " + str(len(solutions)) + " solutions!")
print_intermediate_states(solutions)
```

Now we have selected all amplitudes that involve **f** states. At this point we are all set to generate some data using this amplitude model!


```python
xml_generator = HelicityDecayAmplitudeGeneratorXML()
xml_generator.generate(solutions)
xml_generator.write_to_file('model.xml')
```

#### 5. Creating a data sample


```python
# pycompwa is the python interface to ComPWA's c++ modules
import pycompwa as pwa

# Create particle list
particle_list = pwa.PartList()
with open('model.xml', 'r') as content_file:
    model_file_contents = content_file.read()
    pwa.read_particles(particle_list, model_file_contents)

# Create kinematics
kin = pwa.HelicityKinematics(particle_list, 'model.xml')

# Generate phase space sample
gen = pwa.RootGenerator(particle_list, kin, 12345)
phspSample = pwa.generate_phsp(100000, gen)

# Create Amplitude
with open('model.xml', 'r') as content_file:
    model_file_contents = content_file.read()
    intensity = pwa.incoherent_intensity(model_file_contents, 
                                         particle_list,
                                         kin, phspSample,
                                         phspSample)

# Generate Data
sample = pwa.generate(5000, kin, gen, intensity)
```

#### 6. Plotting
Let's go ahead and make a Dalitz plot of the generated data. First we create a ROOT file containing all of the information inside a TTree.


```python
kin.create_all_subsystems()
rootpl = pwa.RootPlotData(kin, intensity)
rootpl.set_data(sample)
rootpl.write("tree", "rootplot.root", "RECREATE")
```

ComPWA ships with a little plotting module to help you read in ROOT TTree's and generate some common plots using matplotlib


```python
# Plotting
from Plotting.plot import (
    make_dalitz_plots
)
from Plotting.ROOT.rootplotdatareader import open_compwa_plot_data

plot_data = open_compwa_plot_data("rootplot.root")

data_variables = list(plot_data.data.dtype.names)
print("found data variables:", data_variables)
```


```python
#binned_dists = make_binned_distributions(plot_data, var_names)
make_dalitz_plots(plot_data, ['mSq_3_4_vs_2', 'mSq_2_4_vs_3'], bins=50)
```

#### 7. Fitting

All parameter defined and used by the **Intensity** object, can be obtain for it by using the `parameters()` function. Just pass it an empty `ParameterList` object.


```python
par_list = pwa.ParameterList()
intensity.parameters(par_list)
fit_parameters = par_list.get_fit_parameters()
```

Let's save the true parameters in a dictionary so we can compare the fitted values later on. Notice that the `get_fit_parameters()` returns a special object that behave similar to a python list. The contents of the list are FitParameter objects, that have attributes `name, value, error, is_fixed`. The name and error attributes are read only.


```python
true_parameters = {x.name: x.value for x in fit_parameters if not x.is_fixed}
print(true_parameters)
```

To make the fit a bit more interesting, we modify one of the parameters to a different initial value then the true value.


```python
idx = fit_parameters.index("Magnitude_J/psi_to_f2(1270)_0+gamma_-1;f2(1270)_to_pi0_0;")
print("before:", fit_parameters[idx])
fit_parameters[idx].value = 2.0
print("after:",fit_parameters[idx])
# we can also fix or free parameters here
fit_parameters[fit_parameters.index(
    'Phase_J/psi_to_f2(1270)_0+gamma_-1;f2(1270)_to_pi0_0;')].is_fixed = True
print("should be fixed now.... ",fit_parameters[fit_parameters.index(
    'Phase_J/psi_to_f2(1270)_0+gamma_-1;f2(1270)_to_pi0_0;')])
```

Now it's time to start up a set up a fit, which is quite simply.
1. First create an estimator instance of your choice, here a minimum log likelihood (`MinLogLH`). Notice that we use the function tree feature. This create a full evaluation tree, caching the data and the intensity. It can greatly enhance the fit performance!
2. Then create an optimizer instance of your choice, here Minuit2 (`MinuitIF`).


```python
esti = pwa.MinLogLH(kin, intensity, sample, phspSample, phspSample)
esti.enable_function_tree(True)
esti.log_function_tree()

minuitif = pwa.MinuitIF(esti, par_list)
minuitif.enable_hesse(True)

result = minuitif.minimize(par_list)
```

Let's check if the fit parameters are "close to" the true values


```python
fitresult_parameters = {x.name: (x.value, x.error) for x in fit_parameters if not x.is_fixed}
for key, value in fitresult_parameters.items():
    print(key, " fit result:", "{0:.3f}".format(value[0]), "+-", 
          "({0:.3f},".format(value[1][0]), "{0:.3f})".format(value[1][1]),
          " true:", "{0:.3f}".format(true_parameters[key])
         )
```

That's it. You can check some of the other examples to learn about more detailed features of ComPWA.

And give us feedback or contribute ;)!

## Installation
A detailed guide can be found [here](https://github.com/ComPWA/ComPWA/wiki/Installation). The installation basically boils down to:
```bash
git clone --depth 1 https://github.com/ComPWA/ComPWA.git
cd ComPWA
git submodule init
git submodule update
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=<COMPWA_INSTALL_PATH> .. # -DCMAKE_INSTALL_PREFIX is optional
make -j2
make install # this step is optional
```

## Documentation
The [complete documentation](https://compwa.readthedocs.io/en/latest) of the current master branch and specific tags versions are built with readthedocs.
The master branch is automatically built using TravisCI. Probably it is interesting to check out the [log file](https://travis-ci.org/ComPWA/ComPWA) and the projects TravisCI configuration file [<code>.travisCI.yml</code>](https://github.com/ComPWA/ComPWA/blob/master/.travis.yml).
