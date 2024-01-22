
import os
import shutil
import sys
import h5py
import time
from typing import Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
from tqdm import tqdm


ALL_STREAM_OBSERVABLES = ('phi1', 'phi2', 'dist', 'pm1', 'pm2', 'vr')


def read_observables(dir_path, num_max_file=1000):
    """ Read all observables in a directory """

    data = {k: [] for k in ALL_STREAM_OBSERVABLES}

    for i in range(num_max_file):
        data_path = os.path.join(dir_path, f"final_observables/{i}.hdf5")
        if not os.path.exists(data_path):
            continue
        with h5py.File(data_path, 'r') as f:
            for key in data.keys():
                data[key].append(f[key][()])

    # read in the labels
    labels = pd.read_csv(os.path.join(dir_path, 'labels.csv'))

    return data, labels


def write_dataset(path, data):
    """ Save dataset to HDF5 """
    # get the length of each sample and create a pointer
    sample_len = np.array([len(d) for d in data['phi1']])
    ptr = np.cumsum(sample_len)
    ptr = np.insert(ptr, 0, 0)

    # create the dataset
    with h5py.File(path, 'w') as f:
        for key in data.keys():
            f.create_dataset(
                key, data=np.concatenate(data[key]), compression='gzip')
        f.create_dataset('ptr', data=ptr, compression='gzip')


if __name__ == "__main__":
    """ Combine final observables into a single dataset """

    root = "/mnt/ceph/users/tnguyen/stream_sbi/raw_datasets/2params"
    output_name = os.path.join(root, "2params.hdf5")
    output_label_name = os.path.join(root, "2params_labels.csv")
    num_max_dir = 100
    num_max_file = 1000  # number of files per directory

    data = {k: [] for k in ALL_STREAM_OBSERVABLES}
    labels = []

    for i in tqdm(range(num_max_dir), desc='Combining datasets'):
        dir_path = os.path.join(root, f"{i}")
        if not os.path.exists(dir_path):
            continue
        d, t = read_observables(dir_path, num_max_file=num_max_file)
        for key in data.keys():
            data[key] += d[key]
        labels.append(t)

    # concatenate labels
    labels = pd.concat(labels, ignore_index=True)
    labels['pid'] = labels.index

    # write the dataset and labels
    print('Writing dataset and labels...')
    write_dataset(output_name, data)
    labels.to_csv(output_label_name, sep=',', index=False)
