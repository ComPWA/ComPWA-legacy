Quickstart
==========

If you want to use the python interface to ComPWA and/or use the python
modules of ComPWA, setting up a virtual environment (venv) is highly 
recommended. Below are the setup instructions.

Installation & Setup
--------------------

.. note::
   
   Make sure the following software is installed on your system:
   
   * git (optional, for easier updates and if you want to contribute)
   * cmake ( > 3.3 )
   * gcc (> 5.1) or clang
   * `Boost <http://www.boost.org/users/download/>`_\ , version >= 1.54

   Not required, but recommended:

   * python3 + virtualenv (for ComPWA expert system and python interface as well as a python plotting module)
   * `ROOT <http://root.cern.ch/drupal/content/downloading-root>`_\ , version 5.34, 6.08

To install, simply run

.. code-block:: shell

   * `git clone https://github.com/ComPWA/ComPWA.git <COMPWA_SOURCE_PATH>`
   * `cd <COMPWA_SOURCE_PATH> && git submodule init && git submodule update`
   * `mkdir build && cd build`
   * `cmake ../<COMPWA_SOURCE_PATH>`
   * `make`

.. note::

   Here `<COMPWA_SOURCE_PATH>` points to the ComPWA source directory.

**Setup a python virtual environment**

.. code-block:: shell

   virtualenv -p python3 <PATH_OF_YOUR_VENV>
   source <PATH_OF_YOUR_VENV>/bin/activate
   pip install virtualenvwrapper
   source virtualenvwrapper.sh
   add2virtualenv <COMPWA_SOURCE_PATH>/Physics/ExpertSystem
   add2virtualenv <COMPWA_SOURCE_PATH>/Tools
   add2virtualenv <COMPWA_BUILD_DIR>/Tools/PythonInterface

.. note::
   Replace **$PATH_OF_YOUR_VENV** with the path where the venv should be installed.

   `<COMPWA_BUILD_DIR>` points to the ComPWA build directory.
   

**Install requirements for modules**
  
Each python module of ComPWA contains a requirements.txt file. If you want to
use this module simply install the requirements by executing:
  
.. code-block:: shell

   pip install -r <PATH_TO_COMPWA_PYTHON_MODULE>/requirements.txt
    
For example: ``pip install -r Physics/ExpertSystem/requirements.txt``
(assuming you are in the `<COMPWA_SOURCE_PATH>` directory)


Running
-------

This code follows :ref:`this Example <examples-jpsi-to-gammapi0pi0>`,
but was slimlined to be as short as possible. It models the 
:math:`J/\Psi \rightarrow \gamma \pi^0 \pi^0` decay.

1. Creating a model
"""""""""""""""""""

Intensity/Amplitude models can easily be created using the 
:ref:`ComPWA Expert System <compwa-expert-system>`.

.. code:: ipython3

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

    initial_state = [("J/psi", [-1, 1])]
    final_state = [("gamma", [-1, 1]), ("pi0", [0]), ("pi0", [0])]
    
    tbd_manager = StateTransitionManager(initial_state, final_state,
                                         formalism_type='helicity',
                                         topology_building='isobar')
    tbd_manager.set_allowed_interaction_types(
        [InteractionTypes.Strong, InteractionTypes.EM])

    graph_interaction_settings_groups = tbd_manager.prepare_graphs()

    # particles are found by name comparison; so i.e. f2 will find all f2's and f all f's
    tbd_manager.allowed_intermediate_particles = ['f']
    
    (solutions, violated_rules) = tbd_manager.find_solutions(
            graph_interaction_settings_groups)
    
    print("found " + str(len(solutions)) + " solutions!")
    print_intermediate_states(solutions)

.. note::
   The ``StateTransitionManager`` (STM) is the main user interface class of the
   ComPWA expert system. The boundary conditions of your physics problem are 
   defined here, such as the initial state, final state, formalism type, ...

   * ``prepare_graphs()`` of the STM creates all topology graphs, here using
     the isobar model (two-body decays). Also it initializes the graphs with
     the initial and final state and the a set of conservation laws at each
     interaction node.

   * By default all three (strong, EM, weak) interaction types are used in the
     preparation stage. However it is also possible to globally choose the
     allowed interaction types via ``set_allowed_interaction_types()``.

   After the preparation step, you can modifiy the settings returned by
   ``prepare_graphs()`` to your liking. Since this output is quite a lot of
   information, the expertsystem ui is supposed to aid in the configuration
   (especially the STM).
   
   * A subset of particles that are allow as intermediate states can also be
     specified in the STM. Either in the ``init()`` of the STM or setting the
     instance attribute ``allowed_intermediate_particles``.

Now we have created an intensity using the helicity formalism, which includes
several **f** states. At this point we are all set to generate some data using
this amplitude model!

Export the amplitude to a xml file.

.. code:: ipython3

    xml_generator = HelicityDecayAmplitudeGeneratorXML()
    xml_generator.generate(solutions)
    xml_generator.write_to_file('model.xml')

2. Creating data samples
""""""""""""""""""""""""

.. code:: ipython3

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

.. note::
   ``pycompwa`` is the python interface to ComPWA's c++ modules. Read more
   about this :ref:`here <python-ui>`.

   Three important pieces for evaluating an intensity are:

   * The **intensity** itself. It was generated previously and stored within
     the xml model file.

   * A **kinematics** instance. It handles the calculation of the kinematic
     variables that are required for the evaluation of the intensity!
     For example in the helicity formalism: :math:`(s,\theta,\phi)`.
   
   * **Data samples**. For mere visualization of the intensity a phase space
     sample is sufficient. It is mandatory for the normalization of the
     intensity. However when performing fits an additional data sample, to
     which the intensity will be compared to, has to be specified.

3. Plotting
"""""""""""

Let's go ahead and make a Dalitz plot of the generated data. Here there are two
possible ways to achieve this.

1. All of the information can be saved into a file (currenty only ROOT files 
   are supported). Afterwards the file can be read in and processed by the
   ComPWA ``Plotting`` python module.

   .. code:: ipython3

      kin.create_all_subsystems()
      rootpl = pwa.RootPlotData(kin, intensity)
      rootpl.set_data(sample)
      rootpl.write("tree", "rootplot.root", "RECREATE")

      # Plotting
      from Plotting.plot import (
          make_dalitz_plots
      )
      from Plotting.ROOT.rootplotdatareader import open_compwa_plot_data

      plot_data = open_compwa_plot_data("rootplot.root")

      data_variables = list(plot_data.data.dtype.names)
      print("found data variables:", data_variables)

      #binned_dists = make_binned_distributions(plot_data, var_names)
      make_dalitz_plots(plot_data, ['mSq_3_4_vs_2', 'mSq_2_4_vs_3'], bins=50)

2. A very convenient way to access the data inside your python code is through
   the direct data interface ``DataPoints()`` of the ``pycompwa`` module. The
   data can also be passed to the ``Plotting`` module.

   .. code:: ipython3

      kin.create_all_subsystems()

      import numpy as np
      # use the direct data point access via DataPoints
      data_points = pwa.DataPoints(phsp_sample_importance, kin)
      data_array = np.array(data_points)
      variable_names = data_points.get_variable_names()
      print(variable_names)

      # Plotting
      from Plotting.plot import (
          make_dalitz_plots, PlotData
      )

      plotdata = PlotData(column_names=variable_names, data_array=data_array)
      plotdata.particle_id_to_name_mapping = data_points.get_finalstate_id_to_name_mapping()
      # plot a 2d histogram
      make_dalitz_plots(plotdata, ["mSq_3_4_vs_2", "mSq_2_4_vs_3"], bins=50)

.. tip::
   
   ComPWA ships with a little python plotting module (``Plotting``) to help you
   create plots, which are often used (angular distributions, Dalitz plots). It
   uses matplotlib as a backend. You can either hand it data files, or feed it
   directly with data.

   Use it! Instead of creating your own visualization scripts.

4. Fitting
""""""""""

All parameters defined and used by the **Intensity** object, can be obtained by
using the ``parameters()`` function. Just pass it an empty ``ParameterList``
object.

.. code:: ipython3

    par_list = pwa.ParameterList()
    intensity.parameters(par_list)
    fit_parameters = par_list.get_fit_parameters()

Let's save the true parameters in a dictionary so we can compare the fitted
values later on. Notice that the ``get_fit_parameters()`` returns a special 
object that behave similar to a python list. The contents of the list are 
``FitParameter`` objects, with the attributes ``name, value, error, is_fixed``.
The name and error attributes are read only.

.. code:: ipython3

    true_parameters = {x.name: x.value for x in fit_parameters if not x.is_fixed}
    print(true_parameters)

To make the fit a bit more interesting, we modify one of the parameters to a
different initial value then the true value.

.. code:: ipython3

    idx = fit_parameters.index("Magnitude_J/psi_to_f2(1270)_0+gamma_-1;f2(1270)_to_pi0_0;")
    print("before:", fit_parameters[idx])
    fit_parameters[idx].value = 2.0
    print("after:",fit_parameters[idx])
    # we can also fix or free parameters here
    fit_parameters[fit_parameters.index(
        'Phase_J/psi_to_f2(1270)_0+gamma_-1;f2(1270)_to_pi0_0;')].is_fixed = True
    print("should be fixed now.... ",fit_parameters[fit_parameters.index(
        'Phase_J/psi_to_f2(1270)_0+gamma_-1;f2(1270)_to_pi0_0;')])

Now it's time to start up a set up a fit, which is quite simple.

1. First create an estimator instance of your choice, here a minimum log
   likelihood (``MinLogLH``). Notice that we use the function tree feature.
   This create a full evaluation tree, caching the data and the intensity. It
   can greatly enhance the fit performance!
2. Then create an optimizer instance of your choice, here Minuit2
   (``MinuitIF``).

.. code:: ipython3

    esti = pwa.MinLogLH(kin, intensity, sample, phspSample, phspSample)
    esti.enable_function_tree(True)
    esti.log_function_tree()
    
    minuitif = pwa.MinuitIF(esti, par_list)
    minuitif.enable_hesse(True)
    
    result = minuitif.minimize(par_list)

Let's check if the fit parameters are "close to" the true values

.. code:: ipython3

    fitresult_parameters = {x.name: (x.value, x.error) for x in fit_parameters if not x.is_fixed}
    for key, value in fitresult_parameters.items():
        print(key, " fit result:", "{0:.3f}".format(value[0]), "+-", 
              "({0:.3f},".format(value[1][0]), "{0:.3f})".format(value[1][1]),
              " true:", "{0:.3f}".format(true_parameters[key])
             )

This should get you started. You can check some of the other examples to learn
about more detailed features of ComPWA.

And we would be happy to recieve some feedback or contributions ;)!