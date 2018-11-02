def open_compwa_plot_data(input_file_path):
    from ROOT import gROOT, TFile, gDirectory
    from Plotting.plot import PlotData
    pd = PlotData()

    gROOT.Reset()
    # open file
    file = TFile(input_file_path, 'READ')
    file.cd("final_state_id_to_name_mapping")
    for k in gDirectory.GetListOfKeys():
        pd.particle_id_to_name_mapping[k.ReadObj(
        ).GetVal()] = k.ReadObj().GetName()

    from root_numpy import list_trees
    trees = list_trees(input_file_path)
    if "tree_data" in trees:
        pd.data = load_ttree(input_file_path, "tree_data")
    if "tree_weighted_phsp_MC" in trees:
        pd.fit_result_data = load_ttree(
            input_file_path, "tree_weighted_phsp_MC")

    return pd


def load_ttree(fname, tree=None, patterns=None, *kargs, **kwargs):
    """
    Loads a root ttree into a numpy record array
    Further *kargs and *kwargs are passed to root_numpy's root2array.
    >>> df = load_ttree('test.root', 'MyTree', patterns=['x_*', 'y_*'], 
                         selection='x_1 > 100')
    """
    from root_numpy import root2array, list_trees
    import numpy as np

    if tree is None:
        branches = list_trees(fname)
        if len(branches) == 1:
            tree = branches[0]
        else:
            raise ValueError('More than one tree found in {}'.format(fname))

    if not patterns:
        all_vars = None
    else:
        # TODO: select only specifiy branches based on the patterns
        all_vars = None

    arr = root2array(fname, tree, all_vars, *kargs, **kwargs)
    return np.rec.array(arr)
