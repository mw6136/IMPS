import numpy as np
import h5py

class Dataset:

    def __init__(self):
        self.fields = None
        self.dim = None
        self.data = None

def load(fname):

    ds = Dataset()

    dfile = h5py.File(fname, "r")

    # Set the available fields of the dataset
    fields = []

    for key in dfile.keys():
        fields.append(key)

    ds.fields = fields

    # Set the dim of the dataset
    ds.dim = tuple(dfile["Dim"])

    data = {}
    for field in ds.fields:
        if field not in ["Dim", "Info", "TimeInfo"]:
            data[field] = np.array(dfile[field]).reshape((ds.dim[0],ds.dim[1]))[1:-1,1:-1]
        else:
            data[field] = np.array(dfile[field])

    ds.data = data

    return ds