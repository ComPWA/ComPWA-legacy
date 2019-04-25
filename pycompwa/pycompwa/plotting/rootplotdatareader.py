import uproot


def open_compwa_plot_data(input_file_path):
    from pycompwa.plotting import PlotData
    pd = PlotData()

    # open file
    file = uproot.open(input_file_path)
    trees = file.keys()

    file = file.get("final_state_id_to_name_mapping")
    for k, v in file.items():
        pd.particle_id_to_name_mapping[v] = k.decode()[:k.decode().find(';')]

    if "data" in [x.decode()[:x.decode().find(';')] for x in trees]:
        pd.data = load_ttree(input_file_path, "data")
    if "intensity_weighted_phspdata" in [x.decode()[:x.decode().find(';')]
                                         for x in trees]:
        pd.fit_result_data = load_ttree(
            input_file_path, "intensity_weighted_phspdata")

    return pd


def load_ttree(filename, treename, branchnames=None):
    """
    Loads a root ttree into a numpy record array
    If branchnames is None all branches are read.
    """
    import numpy as np
    tree = uproot.open(filename)[treename]

    if not branchnames:
        branchnames = tree.keys()
    array_dict = tree.arrays(branchnames)
    return np.rec.fromarrays(array_dict.values(),
                             names=[x.decode() for x in array_dict.keys()])
