import numpy as np
import scipy
import matplotlib.pyplot as plt
import glob
import os
from typing import Tuple, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor

from sbi import utils as utils
from sbi import analysis as analysis
from sbi.inference.base import infer
from sbi import inference
from sbi.inference import SNPE, simulate_for_sbi, prepare_for_sbi
import pickle

from  subhalo_impact import chi_eval, initial_stream
from rotation_matrix import obs_from_pos6d
from importlib import reload
reload(initial_stream)
import rnn
reload(rnn)
import temp


def filter_nans(stats, empty_bin_indices):
    stats = np.copy(stats)
    new_stats = np.delete(stats, empty_bin_indices)
    return new_stats

def replace_with_neighboring_value(values, values_stats, empty_bin_indices, phi1, phi1_edges):
    values_stats = np.copy(values_stats)
    for i in empty_bin_indices:
        neighboring_bin1 = i-1
        neighboring_bin2 = i+1
        while neighboring_bin1 in empty_bin_indices:
            neighboring_bin1 -= 1
        
        while neighboring_bin2 in empty_bin_indices:
            neighboring_bin2 += 1
        
        median, _, _ = scipy.stats.binned_statistic(phi1, values, statistic='median', bins = 1, range = (phi1_edges[neighboring_bin1], phi1_edges[neighboring_bin2+1]))
        values_stats[i] = median[0]
    return values_stats

def replace_with_neighboring_value_std(std_stats, empty_bin_indices):
    std_stats = np.copy(std_stats)
    for i in empty_bin_indices:
        neighboring_bin1=i-1
        neighboring_bin2=i+1
        while neighboring_bin1 in empty_bin_indices:
            neighboring_bin1 -= 1
        while neighboring_bin2 in empty_bin_indices:
            neighboring_bin2 += 1
        std_stats[i] = np.median([std_stats[neighboring_bin1], std_stats[neighboring_bin2]])
    return std_stats

def calculate_summary_statistics(phi1,phi2,dist,pm1,pm2,vr, bins):
    # bins = np.arange(-25, 25, 0.5)
    # bins = 100
    phi2_median, phi1_edges, _ = scipy.stats.binned_statistic(phi1, phi2, statistic='median', bins = bins)
    phi2_std, _, _ = scipy.stats.binned_statistic(phi1, phi2, statistic='std', bins = bins)
    count, _, _ = scipy.stats.binned_statistic(phi1, phi2, statistic='count', bins = bins)    
    vr_median, _, _ = scipy.stats.binned_statistic(phi1, vr, statistic='median', bins = bins)
    vr_std, _, _ = scipy.stats.binned_statistic(phi1, vr, statistic='std', bins = bins)
    pm1_median, _, _= scipy.stats.binned_statistic(phi1, pm1, statistic='median', bins = bins)
    pm1_std, _, _ = scipy.stats.binned_statistic(phi1, pm1, statistic='std', bins = bins)
    pm2_median, _, _= scipy.stats.binned_statistic(phi1, pm2, statistic='median', bins = bins)
    pm2_std, _, _ = scipy.stats.binned_statistic(phi1, pm2, statistic='std', bins = bins)
    dist_median, _, _ = scipy.stats.binned_statistic(phi1, dist, statistic='median', bins = bins)
    dist_std, _, _ = scipy.stats.binned_statistic(phi1, dist, statistic='std', bins = bins)

    # '''fill in empty bins with neighboring values'''
    # nan_indices = np.argwhere(np.isnan(phi2_median)).flatten()

    # count_new = count
    # phi1_edges_new = phi1_edges
    # phi2_median_new = replace_with_neighboring_value(phi2, phi2_median, nan_indices, phi1, phi1_edges)
    # vr_median_new = replace_with_neighboring_value(vr, vr_median, nan_indices, phi1, phi1_edges)
    # pm1_median_new = replace_with_neighboring_value(pm1, pm1_median, nan_indices, phi1, phi1_edges)
    # pm2_median_new = replace_with_neighboring_value(pm2, pm2_median, nan_indices, phi1, phi1_edges)
    # dist_median_new = replace_with_neighboring_value(dist, dist_median, nan_indices, phi1, phi1_edges)

    # phi2_std_new = replace_with_neighboring_value_std(phi2_std, nan_indices)
    # vr_std_new = replace_with_neighboring_value_std(vr_std, nan_indices)
    # pm1_std_new = replace_with_neighboring_value_std(pm1_std, nan_indices)
    # pm2_std_new = replace_with_neighboring_value_std(pm2_std, nan_indices)
    # dist_std_new = replace_with_neighboring_value_std(dist_std, nan_indices)

    # '''filter out nan values, resulting in different sizes of output'''
    # nan_indices = np.argwhere(np.isnan(phi2_median)).flatten()
    # phi1_edges_new = filter_nans(phi1_edges, nan_indices)
    # count_new = filter_nans(count, nan_indices)
    # phi2_median_new = filter_nans(phi2_median, nan_indices)
    # vr_median_new = filter_nans(vr_median, nan_indices)
    # pm1_median_new = filter_nans(pm1_median, nan_indices)
    # pm2_median_new = filter_nans(pm2_median, nan_indices)
    # dist_median_new = filter_nans(dist_median, nan_indices)
    # phi2_std_new = filter_nans(phi2_std, nan_indices)
    # vr_std_new = filter_nans(vr_std, nan_indices)
    # pm1_std_new = filter_nans(pm1_std, nan_indices)
    # pm2_std_new = filter_nans(pm2_std, nan_indices)
    # dist_std_new = filter_nans(dist_std, nan_indices)

    # output = np.vstack((phi1_edges_new[:-1], phi2_median_new, count_new/len(phi1), 
    #                     vr_median_new, pm1_median_new, pm2_median_new, 
    #                     dist_median_new, phi2_std_new, vr_std_new, pm1_std_new, 
    #                     pm2_std_new, dist_std_new))
    
    output = np.vstack((phi1_edges[:-1], phi2_median, count/len(phi1), 
                        vr_median, pm1_median, pm2_median, 
                        dist_median, phi2_std, vr_std, pm1_std, 
                        pm2_std, dist_std))
    output = torch.from_numpy(output)
    return output

def filter_empty_sim (pid_list, bins): 
    new_list = list()
    for pid in pid_list:
        parameters, phi1, phi2, dist, pm1, pm2, vr = temp.read_observables_hdf5('sim'+str(pid)+'.hdf5')
        phi2_median, phi1_edges, _ = scipy.stats.binned_statistic(phi1, phi2, statistic='median', bins = bins)
        nan_indices = np.argwhere(np.isnan(phi2_median)).flatten()
        if len(nan_indices) == 0 : 
            new_list.append(pid)
    return new_list

def find_sims_pid (pid_range, M_limit, bins):
    '''
    pid_range : tuple of ints (start, end) for varying different parameters
    M_limit : float, max M_sat value
    bins : np array for phi1 ranges
    '''
    start, end = pid_range
    pid_list = list()
    with open("pid.txt") as f:
        text = f.read()
    rows = text.split('\n')
    rows = rows[(start + 1) : (end + 1)]
    for r in rows :
        data = r.split(",")
        m = data[5]
        if float(m) < M_limit :
            pid_list.append(int(data[0]))
    pid_list = filter_empty_sim(pid_list, bins)
    return np.array(pid_list)


def simulator (params, scatter = 0.1, seqlen = 100):
    'simulation output saved in sim0.hdf5'
    params = np.array(params, dtype=float) 
    log_Msat = params[0]
    M_sat = 10**log_Msat
    vz = params[1]

    pid = 0
    r = 0.2
    phi = 250 
    vphi = 35  
    t_a = 0.2 
    phi_a = -4 
    rs_sat = 1.05 * (M_sat*10*10)**0.5 
    tmax=4 

    chi_eval(r,phi,vphi,vz,M_sat,tmax,t_a,phi_a,rs_sat,pid)
    temp.save_observables_hdf5(r,phi,vphi,vz,M_sat,tmax,t_a,phi_a,rs_sat,pid)
    parameters, phi1,phi2,dist,pm1,pm2,vr = temp.read_observables_hdf5(f'sim{pid}.hdf5')
    output = calculate_summary_statistics(phi1,phi2,dist,pm1,pm2,vr)
    output = (output + torch.randn_like(output) * scatter).T
    output = F.pad(output, (0, 0, 0, seqlen-len(output)), value = 0.0)
    return output

def load_non_uniform(pid_list, s_len, bins):
    'pid number list to file name list'
    pid_list = pid_list.astype(str)
    pid_list = np.char.add('sim',pid_list)
    pid_list = np.char.add(pid_list, '.hdf5')

    my_files = pid_list
    batch_size = len(my_files)
    # seq_len = torch.zeros(batch_size)

    # input: num_sim * 2
    # outpout: num_sim * 12
    theta = torch.zeros(batch_size, 2)
    's_len + 1: additional space to store seqlen'
    x = torch.zeros(batch_size, s_len, 12)

    i = 0
    '''load the first? num_sim files from observables'''
    for hf in my_files:
        parameters,phi1,phi2,dist,pm1,pm2,vr = temp.read_observables_hdf5(hf)
        output = calculate_summary_statistics(phi1,phi2,dist,pm1,pm2,vr,bins)
        output = (output + torch.randn_like(output) * 0.1).T
        len_output = len(output)
        if len_output != 0:
            # seq_len[i] = len_output
            pad_len = s_len - len_output
            # print(hf, ':',pad_len, "empty?")
            '''need to truncate output if len > seq_len?'''
            output = F.pad(output, (0, 0, 0, pad_len), value = 0.0)

            # sequence = torch.randn(1)
            # sequence[0] = seq_len[i]
            # sequence = sequence.repeat(s_len, 1)
            'append original output to it'
            # output = torch.cat([sequence, output], dim=1)
            x[i] = output
            # if i == 0:
                # print(output)

            # log_Msat, vz from dataset "parameters"
            log_Msat = np.log10(parameters[4])
            vz = parameters[3]
            theta[i] = torch.as_tensor([vz, log_Msat])

            i += 1
        # if i >= num_sim: 
            # print(x.shape, "load sims")
            # seq_len = np.array(seq_len)
            # print('0 seqlen', np.all(seq_len))
            # break
    return theta, x

if __name__ == '__main__':
    num_dim = 2
    prior_min = [-50, -5]
    prior_max = [0, -1]
    prior_uniform = utils.BoxUniform(low = torch.as_tensor(prior_min), 
                                    high = torch.as_tensor(prior_max))
    prior_normal = torch.distributions.MultivariateNormal(5 * torch.ones(num_dim), 
                                                        1 * torch.eye(num_dim))
    # simulator_wrapper, prior = prepare_for_sbi(simulator, prior_uniform)

    n_in, n_out, num_layers, hidden_features, seqlen = 12, 2, 2, 8, 100
    embedding_net = rnn.EmbeddingNet(n_in, n_out, hidden_features, num_layers)
    neural_posterior = utils.posterior_nn(
        model="maf", embedding_net=embedding_net, hidden_features=8)
    inference = SNPE(prior=prior_uniform, density_estimator=neural_posterior)

    bins = np.arange(-25, 25, 0.5)
    pid_range1 = (0, 1700)
    pid_range2 = (2200, 4370)
    M_limit = 0.1 # in unit of M_sun
    pid_list1 = find_sims_pid (pid_range1, M_limit, bins)
    pid_list2 = find_sims_pid (pid_range2, M_limit, bins)
    pid_list = np.hstack((pid_list1, pid_list2))

    theta, x = load_non_uniform(pid_list, seqlen, bins)
    # print('theta', theta)
    # print('x', x)
    density_estimator = inference.append_simulations(theta, x).train()
    posterior = inference.build_posterior(density_estimator)
    with open("my_posterior.pkl", "wb") as handle:
        pickle.dump(posterior, handle)