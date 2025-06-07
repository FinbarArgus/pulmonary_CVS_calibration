import os, sys
root_dir_path = os.path.join(os.path.dirname(__file__), '../')
parent_dir_path = os.path.join(os.path.dirname(__file__), '../../')
sys.path.append(os.path.join(parent_dir_path, 'python_ECG_tools'))
sys.path.append(os.path.join(parent_dir_path, 'ultrasound_uncertainty_tools'))
sys.path.append(os.path.join(parent_dir_path, 'circulatory_autogen'))
import pandas as pd
import numpy as np
import time
import json
from matplotlib import pyplot as plt

from haemotools.BB_echo_functions import BBdata
from haemotools.signal_analysis import signal_tools
from read_tomtec_worksheets import read_tomtec_worksheets
from get_valve_data import write_to_json_file
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )


def generate_intrabeat_periods_json(patient_num, project_dir):

    do_plot = True

    ec = BBdata()
    st = signal_tools()

    data_dir_path = os.path.join(project_dir, 'data/pulmonary/Alfred_ECG_Pre_Post_P1_13')
    plot_dir = os.path.join(data_dir_path, 'plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
        
    fs = 180 
    dt = 1/fs

    
    ecg_path = os.path.join(data_dir_path, f'P{patient_num}_pre_PEA.csv')
            
    # read in ECG data with numpy
    ecg_data = np.genfromtxt(ecg_path, delimiter=',', skip_header=1)
    time = np.arange(0, len(ecg_data)*dt, dt)

    if do_plot:
        # plot ecg data
        fig, ax = plt.subplots(1,1)
        plt.plot(time, ecg_data)
        plt.savefig(os.path.join(plot_dir, 'ecg_data.png'))

    time_segs, ecg_segs, quality, heart_period_dict =  \
            st.slice_by_ecg_with_quality(time, ecg_data, fs=fs)
    
    num_segs = len(time_segs)

    if num_segs < 1:
        print(f' data is too short')
        exit()

    # TODO is the below the best way to do it
    # temporarily fix the ventricle relaxation period to 0.2*T
    for II in range(len(heart_period_dict['T'])):
        heart_period_dict['T_vr'][II] = heart_period_dict['T'][II]*0.2
    print(heart_period_dict['T_vr'])


    quality_threshold = 0.5
    print(quality)
    idxs_to_remove_ecg = [II for II in range(num_segs) if (quality[II] < quality_threshold)]
    print('due to quality, removing indices: ', idxs_to_remove_ecg)

    time_segs = np.delete(time_segs, idxs_to_remove_ecg)
    ecg_segs = np.delete(ecg_segs, idxs_to_remove_ecg, axis=0)
    heart_period_dict['T'] = np.delete(heart_period_dict['T'], idxs_to_remove_ecg).tolist()
    heart_period_dict['T_vc'] = np.delete(heart_period_dict['T_vc'], idxs_to_remove_ecg).tolist()
    heart_period_dict['T_vr'] = np.delete(heart_period_dict['T_vr'], idxs_to_remove_ecg).tolist()
    heart_period_dict['T_ac'] = np.delete(heart_period_dict['T_ac'], idxs_to_remove_ecg).tolist()
    heart_period_dict['T_ar'] = np.delete(heart_period_dict['T_ar'], idxs_to_remove_ecg).tolist()
    heart_period_dict['t_astart'] = np.delete(heart_period_dict['t_astart'], idxs_to_remove_ecg).tolist()

    num_segs_reduced = len(time_segs)

    mean_T = np.mean(heart_period_dict['T'])
    mean_T_vc = np.mean(heart_period_dict['T_vc'])
    mean_T_ac = np.mean(heart_period_dict['T_ac'])
    mean_T_ar = np.mean(heart_period_dict['T_ar'])
    mean_t_astart = np.mean(heart_period_dict['t_astart'])

    opt_idx = np.argmin([(heart_period_dict['T'][II]-mean_T)**2 +
                    (heart_period_dict['T_vc'][II]-mean_T_vc)**2 +
                    (heart_period_dict['T_ac'][II]-mean_T_ac)**2 +
                    (heart_period_dict['T_ar'][II]-mean_T_ar)**2 +
                    (heart_period_dict['t_astart'][II]-mean_t_astart)**2
                    for II in range(num_segs_reduced)])

    opt_T = heart_period_dict['T'][opt_idx]
    opt_T_vc = heart_period_dict['T_vc'][opt_idx]
    opt_T_ac = heart_period_dict['T_ac'][opt_idx]
    opt_T_ar = heart_period_dict['T_ar'][opt_idx]
    opt_T_vr = heart_period_dict['T_vr'][opt_idx]
    opt_t_astart = heart_period_dict['t_astart'][opt_idx]

    print(f'optimal T: {opt_T}')    
    print(f'optimal T_vc: {opt_T_vc}')
    print(f'optimal T_ac: {opt_T_ac}')
    print(f'optimal T_ar: {opt_T_ar}')
    print(f'optimal T_vr: {opt_T_vr}')   
    print(f'optimal t_astart: {opt_t_astart}')
    optimal_dict = {'T': opt_T, 'T_vc': opt_T_vc, 'T_ac': opt_T_ac, 
                    'T_ar': opt_T_ar, 'T_vr': opt_T_vr, 't_astart': opt_t_astart}
    
    fig, ax = plt.subplots()
    [plt.plot(time_seg,seg) for time_seg, seg in zip(time_segs, ecg_segs)]
    plt.vlines(time_segs[opt_idx][0], min(ecg_segs[opt_idx]), max(ecg_segs[opt_idx]), 'k', 'dashed')
    plt.vlines(time_segs[opt_idx][-1], min(ecg_segs[opt_idx]), max(ecg_segs[opt_idx]), 'k', 'dashed')

    plt.savefig(os.path.join(plot_dir, f'ecg_segs_{patient_num}.png'))
    plt.close()

    no_conv = 1.0
    
    constants_of_interest = [('T', 'T',
                            'second', no_conv, 'Alfred_database'),
                            ('T_vc', 'T_vc',
                            'second', no_conv, 'Alfred_database'),
                            ('T_vr', 'T_vr',
                            'second', no_conv, 'Alfred_database'),
                            ('T_ac', 'T_ac',
                            'second', no_conv, 'Alfred_database'),
                            ('T_ar', 'T_ar',
                            'second', no_conv, 'Alfred_database'),
                            ('t_astart', 't_astart',
                            'second', no_conv, 'Alfred_database')]

    variables_of_interest = []
    
    reduced_data_dict = {f'{patient_num}': optimal_dict}
    save_dir= os.path.join(data_dir_path, '../ground_truth_for_CA')
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    
    write_to_json_file(reduced_data_dict, variables_of_interest, constants_of_interest,
                    save_dir, save_name='alfred_periods', sample_rate=fs)
    
if __name__ == "__main__":
    
    if len(sys.argv) == 3:
        patient_num = sys.argv[1]
        project_dir = sys.argv[2]
        
    else:
        print("usage:  python get_intrabeat_periods.py patient_num project_dir") 
        exit()

    generate_intrabeat_periods_json(patient_num, project_dir)

    
