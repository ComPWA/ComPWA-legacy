from ROOT import gROOT, TFile, gDirectory
from Plotting.plot import PlotData

def open_compwa_plot_data(input_file_path):
    gROOT.Reset()
    
    #open file
    file = TFile(input_file_path, 'READ')

    data_integral = file.Get('data_integral')
    weighted_phsp_integral = file.Get('weighted_phsp_integral')

    pd = PlotData()
    pd.data_integral = data_integral.GetVal()
    pd.weighted_phsp_integral = weighted_phsp_integral.GetVal()
    pd.data = load_root("plot.root", "tree_data")
    pd.fit_result_data = load_root("plot.root", "tree_weighted_phsp_MC")
    file.cd("final_state_id_to_name_mapping")
    for k in gDirectory.GetListOfKeys():
        pd.particle_id_to_name_mapping[k.ReadObj().GetVal()] = k.ReadObj().GetName()
    return pd


def load_root(fname, tree=None, patterns=None, *kargs, **kwargs):
    """
    Loads a root file into a pandas DataFrame.
    Further *kargs and *kwargs are passed to root_numpy's root2array.
    >>> df = load_root('test.root', 'MyTree', patterns=['x_*', 'y_*'], selection='x_1 > 100')
    If the root file contains a branch called index, it will become the DataFrame's index.
    """
    from pandas import DataFrame
    from root_numpy import root2array, list_trees

    if tree == None:
        branches = list_trees(fname)
        if len(branches) == 1:
            tree = branches[0]
        else:
            raise ValueError('More than one tree found in {}'.format(fname))

    if not patterns:
        all_vars = None
    else:
        # index is always loaded if it exists
        patterns.append('index')
        all_vars = get_matching_variables(fname, tree, patterns)

    arr = root2array(fname, tree, all_vars, *kargs, **kwargs)
    if 'index' in arr.dtype.names:
        df = DataFrame.from_records(arr, index='index')
    else:
        df = DataFrame.from_records(arr)
    return df