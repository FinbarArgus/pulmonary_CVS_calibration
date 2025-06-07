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
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )

def get_needed_data_inv_and_tomtec(inv_data_path, tomtec_data_path,
                                        measurements_needed_inv, measurements_needed_tomtec, 
                                        nice_names_inv=None, nice_names_tomtec=None,
                                        measurements_needed_inv_extra=None):
    """
    Creates a dict from the path to pickled invasive data file and csv tomtec data
    This dict only has the patients with all of the measurements in 
    measurements_needed_inv and measurements_needed_tomtec

    measurements_needed_inv_extra is an array of measurements that is required for each measurement
    in measurements_needed_inv
    
    """

    if nice_names_inv in [[], None]:
        nice_names_inv = measurements_needed_inv
    if nice_names_tomtec in [[], None]:
        nice_names_tomtec= measurements_needed_tomtec
    if measurements_needed_inv_extra in [[], None]:
        mesurements_needed_inv_extra = []
    nice_names_all = nice_names_inv + [measurements_needed_inv_extra[II]+'_'+nice_names_inv[JJ] for 
                II in range(len(measurements_needed_inv_extra)) for JJ in range(len(nice_names_inv))] + nice_names_tomtec
    
    # get inv data dict
    inv_dict = pd.read_pickle(inv_data_path)
    all_labels = ec.dict_counter(inv_dict)

    #get tomtec (echo) dataframe
    df = pd.read_csv(tomtec_data_path)
    df['patient_ID'].replace({'_':''}, regex=True, inplace=True)
    desc_df = pd.read_csv(tomtec_desc_path) 

    # get the patients that have the measurements we need
    patients_with_measurement = []
    for measurement_name in measurements_needed_inv:
        patients_with_measurement.append(all_labels[measurement_name]['idxs'])

    for measurement_name in measurements_needed_tomtec:
        patients_with_measurement.append([entry for entry in 
            df[~np.isnan(df[measurement_name])]['patient_ID'].to_list()])

    all_set = set.intersection(*[set(x) for x in patients_with_measurement])
    patients_with_all = list(all_set)
    # patients_with_all = list(set(patients_with_v_aov) & set(patients_with_D_aov)).sort()

    num_patients = len(patients_with_all)

    inv_needed_data = {}
    for patient_id in patients_with_all:
        inv_needed_data[patient_id] = {}
        for measurement_name, nice_name in zip(measurements_needed_inv, nice_names_inv):
            inv_needed_data[patient_id][nice_name] = []
            inv_needed_data[patient_id][nice_name+'_count'] = 0
            inv_needed_data[patient_id][nice_name+'_type'] = 'array'
            # add extra measurements for each needed measurement name (like 'II')
            for extra_measurement in measurements_needed_inv_extra:
                inv_needed_data[patient_id][extra_measurement+'_'+nice_name] = []
                inv_needed_data[patient_id][extra_measurement+'_'+nice_name+'_count'] = 0
                inv_needed_data[patient_id][extra_measurement+'_'+nice_name+'_type'] = 'array'

            for Xnum_idx in inv_dict[patient_id].keys():
                if measurement_name in inv_dict[patient_id][Xnum_idx].keys():
                    inv_needed_data[patient_id][nice_name].append(inv_dict[patient_id][Xnum_idx][measurement_name])
                    inv_needed_data[patient_id][nice_name+'_count'] += 1
                    for extra_measurement in measurements_needed_inv_extra:
                        try:
                            inv_needed_data[patient_id][extra_measurement + '_' + nice_name].append(inv_dict[patient_id][Xnum_idx][extra_measurement])
                        except:
                            print(f'measurement {extra_measurement} was not found with corresponding'
                                  f'measurement {measurement_name} for patient {patient_id}, {Xnum_idx}.'
                                  f'all extra_measurements should'
                                  f'always be available, exiting')
                            exit()
                        inv_needed_data[patient_id][extra_measurement +'_' + nice_name + '_count'] += 1

    needed_data = inv_needed_data

    tomtec_needed_data = df[df['patient_ID'].isin(patients_with_all)][['patient_ID'] + measurements_needed_tomtec]
    for patient_id in patients_with_all:
        for measurement_name, nice_name in zip(measurements_needed_tomtec, nice_names_tomtec):
            needed_data[patient_id][nice_name] = [tomtec_needed_data[tomtec_needed_data['patient_ID']==patient_id]
                    [measurement_name].iloc[0]]
            needed_data[patient_id][nice_name+'_count'] = 1
            needed_data[patient_id][nice_name+'_type'] = 'float'

    return needed_data, patients_with_all, nice_names_all

def write_to_json_file(data_dict, variables_of_interest, constants_of_interest, save_path, 
                       extra_df=None, sample_rate=100, save_name='lv_estimation'):
    """
    variables_of_interest = [(data_variable_name, model_variable_name,
                               variable_type {min, max, mean, series, or [series, min]},
                               state_or_alg, unit, conversion, weighting_for_importance)]

    constants_of_interest = [('data_variable_name', 'model_variable_name')
                            'unit', convergence, 'description_of_measurement')]

    """
    for participant_id in data_dict.keys():
        entry_list= []
        constant_list= []
        # loop through variables of interest and get vals from tomtec_df
        for measurement in data_dict[participant_id]:
            write_type = None
            for tup in variables_of_interest:
                if tup[0] == measurement:
                    variable_tuple = tup
                    write_type = 'observable'
                    break
                else:
                    write_type == None
            for tup in constants_of_interest:
                if tup[0] == measurement:
                    variable_tuple = tup
                    write_type = 'constant'
                    break
                else:
                    write_type == None

            if write_type == 'observable':
                model_variable_name = variable_tuple[1]
                if isinstance(variable_tuple[2], str):
                    variable_type_list = [variable_tuple[2]]
                elif isinstance(variable_tuple[2], list):
                    variable_type_list = variable_tuple[2]
                else:
                    print('3rd variable of interest tuple entry should be string or list of strings')
                    exit()

                unit = variable_tuple[3]
                conversion = variable_tuple[4]
                weight = variable_tuple[5]
                CV = variable_tuple[6]
                name_for_plotting = variable_tuple[7]

                for variable_type in variable_type_list:
                    if variable_type == 'min':
                        value = np.min(data_dict[participant_id][measurement])
                        constant_or_series = 'constant'
                    elif variable_type == 'mean':
                        value = np.mean(data_dict[participant_id][measurement])
                        constant_or_series = 'constant'
                    elif variable_type == 'max':
                        value = np.max(data_dict[participant_id][measurement])
                        constant_or_series = 'constant'
                    elif variable_type == 'series':
                        value = data_dict[participant_id][measurement]
                        constant_or_series = 'series'
                    else:
                        print(f'variable type of {variable_type} is not implemented')
                        exit()

                    entry = {'variable': model_variable_name,
                             'data_type': constant_or_series,
                             'unit': unit,
                             'name_for_plotting': name_for_plotting,
                             'weight': weight,
                             'std': 0.01*CV*value*conversion,
                             'obs_type': variable_type,
                             'value': value*conversion}
                    if variable_type == 'series':
                        entry['sample_rate'] = sample_rate
                    entry_list.append(entry)
                
                # if variable_type = 'series':
                #     # also add min and max of series
                #     entry[data_type] = 'min' 


            elif write_type == 'constant':
                model_variable_name = variable_tuple[1]
                unit = variable_tuple[2]
                conversion = variable_tuple[3]
                data_reference = variable_tuple[4]
                value = data_dict[participant_id][measurement]

                constant_entry = {'variable_name': model_variable_name,
                         'units': unit,
                         'value': value*conversion,
                         'data_reference': data_reference}
                constant_list.append(constant_entry)

        if len(entry_list) > 0:
            json_df = pd.DataFrame(entry_list)
            if extra_df is not None:
                json_df = pd.concat([json_df, extra_df])
            print(json_df)
            # json_tomtec_df.to_json(output_json_file_name)
            result = json_df.to_json(orient='records')
            parsed = json.loads(result)
            with open(os.path.join(save_path, f'{save_name}_observables_{participant_id}.json'), 'w') as wf:
                json.dump(parsed, wf, indent=2)
        
        if len(constant_list) > 0:
            constants_json_df = pd.DataFrame(constant_list)
            constants_result = constants_json_df.to_json(orient='records')
            constants_parsed = json.loads(constants_result)
            with open(os.path.join(save_path, f'{save_name}_constants_{participant_id}.json'), 'w') as wf:
                json.dump(constants_parsed, wf, indent=2)

def get_pressure_segments_and_intrabeat_periods_from_ecg(data_dict, plot_dir):

    pressure_segs_all = {}
    heart_period_all_dict = {}
    for patient_id in all_patient_ids:
        pressure_segs_all[patient_id] = {}
        heart_period_all_dict[patient_id] = {}
        for measurement in nice_names_inv:
            pressure_segs_all[patient_id][measurement] = {}
            pressure_segs_all[patient_id][measurement]['values'] = []
            pressure_segs_all[patient_id][measurement]['quality'] = []
            pressure_segs_all[patient_id][measurement]['ecg'] = []
            heart_period_all_dict[patient_id][measurement] = {}
            for idx in range(data_dict[patient_id][measurement+'_count']):
                ecg = data_dict[patient_id]['II_' + measurement][idx]
                pressure = data_dict[patient_id][measurement][idx]
                if len(pressure) < 3*sample_rate:
                    print(f'patient {patient_id}, {measurement}, {idx} data is too short')
                    continue

                print(f'patient {patient_id}, {measurement}, {idx}')
                pressure_segs, ecg_segs, quality, heart_period_dict = st.slice_by_ecg_with_quality(pressure, ecg)

                if pressure_segs == []:
                    print(f'patient {patient_id}, {measurement}, {idx} data is too short')
                    continue
                    # continue

                # TODO remove the below
                # temporarily fix the ventricle relaxation period to 0.2*T
                for II in range(len(heart_period_dict['T'])):
                    heart_period_dict['T_vr'][II] = heart_period_dict['T'][II]*0.2

                fig, ax = plt.subplots()
                [plt.plot(seg) for seg in pressure_segs]

                plt.savefig(os.path.join(plot_dir, f'pressure_segs_{patient_id}_{measurement}_{idx}.eps'))
                plt.close()

                # for seg_idx in range(len(ecg_segs)):
                #     fig, ax = plt.subplots()
                #     t = np.linspace(0, heart_period_dict['T'][seg_idx], len(ecg_segs[seg_idx]))
                #     ax.plot(t, ecg_segs[seg_idx])
                #     ax.axvline(x=heart_period_dict['T_vc'][seg_idx], color='b', label='T_vc_end')
                #     ax.axvline(x=(heart_period_dict['T_vc'][seg_idx] + heart_period_dict['T_vr'][seg_idx]),
                #                 color='r', label='T_vr_end')
                #     ax.axvline(x=heart_period_dict['t_astart'][seg_idx], color='m', label='t_astart')
                #     ax.axvline(x=(heart_period_dict['t_astart'][seg_idx] + heart_period_dict['T_ac'][seg_idx] +
                #                    heart_period_dict['T_ar'][seg_idx]), color='g', label='T_ar_end')
                #     ax.legend()
                #     plt.savefig(os.path.join(plot_dir, f'ecg_check_{patient_id}_{measurement}_{idx}_{seg_idx}.pdf'))
                #     plt.close()

                fig0, ax0 = plt.subplots()
                ax1 = ax0.twinx()
                ax0.plot(pressure, 'k')
                ax1.plot(ecg, 'r')
                plt.savefig(os.path.join(plot_dir, f'pressure_{patient_id}_{measurement}_{idx}.eps'))
                plt.close()

                # save all pressure_segs to list
                pressure_segs_all[patient_id][measurement]['values'] += pressure_segs
                pressure_segs_all[patient_id][measurement]['quality'] += quality
                pressure_segs_all[patient_id][measurement]['ecg'] += ecg_segs
                for key in heart_period_dict:
                    if key not in heart_period_all_dict[patient_id][measurement].keys():
                        heart_period_all_dict[patient_id][measurement][key] = []
                    heart_period_all_dict[patient_id][measurement][key] += heart_period_dict[key]
    return pressure_segs_all, heart_period_all_dict


if __name__ == '__main__':
    
    ec = BBdata()
    st = signal_tools()

    inv_data_path = '/eresearch/heart/farg967/Sandboxes/Stephen/Biobeat/IVP_data_test.pickle'
    tomtec_data_path = '/eresearch/heart/farg967/Sandboxes/Finbar/tomtec/tomtec_data.csv'
    tomtec_desc_path = '/eresearch/heart/farg967/Sandboxes/Finbar/tomtec/tomtec_desc.csv'
    save_path_combined = '/eresearch/heart/farg967/Sandboxes/Finbar/combined/'

    # plotting
    plot_dir = os.path.join(root_dir_path, 'images')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    # define 
    measurements_needed_inv = ['LV', 'AO']
    measurements_needed_inv_extra = ['II']
    nice_names_inv = measurements_needed_inv
    # nice_names_inv = ['P_lv', 'P_ao'] # THese correspond to the above
    measurements_needed_tomtec = ['US.CA.LVOT.VMAX:ANTFLOW:DOPPLER', 'US.CA.LVOT.DIAM:BMODE', 
            'US.CA.AV.AVA:CONT_VTI:ANTFLOW:DOPPLER', 
            'US.CA.AO.DIAM_STJ:BMODE', 
            'US.CA.LA.VOLES:BP:BMODE', 
            'US.CA.LA.VOLED:CPA2D', 
            'US.CA.LV.ESV:LVA4D', 
            'US.CA.LV.EDV:LVA4D']
            # 'US.CA.RA.VOLES:SP:BMODE:A4C']

    nice_names_tomtec = measurements_needed_tomtec
    # nice_names_tomtec = ['vel_aov', 'd_lvoft', 'd_aov', 'd_ao'] # THese correspond to the above
    sample_rate = 240
    dt = 1/sample_rate

    print(f'getting data for the patients and measurements specified')

    data_dict, all_patient_ids, nice_measurement_names = get_needed_data_inv_and_tomtec(inv_data_path, tomtec_data_path, 
                                                    measurements_needed_inv, 
                                                    measurements_needed_tomtec, 
                                                    nice_names_inv=nice_names_inv, 
                                                    nice_names_tomtec=nice_names_tomtec,
                                                    measurements_needed_inv_extra=measurements_needed_inv_extra)

    num_patients = len(all_patient_ids)

    print(f'data dict for patients and needed measurements created')
    print(f'running analysis with {num_patients} patients')
    print(f'____________________________________________')
    print(f'slicing pressure segments from ecg')

    pressure_segs_all, heart_period_all_dict = get_pressure_segments_and_intrabeat_periods_from_ecg(data_dict, plot_dir)


    print('averaging each patients pressure segments and throwing out outliers')
    quality_threshold = 0.5
    
    # create data dict with only the needed data
    reduced_data_dict = {}
    reduced_heart_period_dict = {}

    for patient_id in all_patient_ids:
        reduced_data_dict[patient_id] = {}
        reduced_heart_period_dict[patient_id] = {}
        for measurement in nice_measurement_names:
            if measurement in nice_names_inv:
                reduced_heart_period_dict[patient_id][measurement] = {}
                num_segs = len(pressure_segs_all[patient_id][measurement]['values'])
                # get mean max pressure and mean min pressure
                pressure_segs = pressure_segs_all[patient_id][measurement]['values']
                if len(pressure_segs) == 0:
                    print(f'Errors found in patient {patient_id}, ignoring their data')
                    del reduced_data_dict[patient_id]
                    del reduced_heart_period_dict[patient_id]
                    break
                quality = pressure_segs_all[patient_id][measurement]['quality']
                idxs_to_remove_ecg = [II for II in range(num_segs) if (quality[II] < quality_threshold )]

                mean_of_max = np.mean([max(pressure_segs[II]) for II in range(num_segs) if II not in idxs_to_remove_ecg])
                std_of_max = np.std([max(pressure_segs[II]) for II in range(num_segs) if II not in idxs_to_remove_ecg])

                idxs_to_remove_pressure = [II for II in range(num_segs) if (max(pressure_segs[II]) > mean_of_max + 2*std_of_max) or \
                                        (max(pressure_segs[II]) < mean_of_max - 2*std_of_max)]

                mean_of_min = np.mean([min(pressure_segs[II]) for II in range(num_segs)])
                std_of_min = np.std([min(pressure_segs[II]) for II in range(num_segs)])

                idxs_to_remove_pressure += [II for II in range(num_segs) if min(pressure_segs[II]) > mean_of_min + 2* std_of_min or \
                                         min(pressure_segs[II]) < mean_of_min - 2* std_of_min]
                
                idxs_to_remove = idxs_to_remove_ecg + idxs_to_remove_pressure
                
                print(f'removing {len(idxs_to_remove)}/{num_segs} for patient {patient_id}, '
                      f'measurement {measurement}')

                # remove
                for key in pressure_segs_all[patient_id][measurement].keys():
                    pressure_segs_all[patient_id][measurement][key] = \
                        np.delete(pressure_segs_all[patient_id][measurement][key], idxs_to_remove)
                for key in heart_period_all_dict[patient_id][measurement].keys():
                    heart_period_all_dict[patient_id][measurement][key] = \
                        np.delete(heart_period_all_dict[patient_id][measurement][key], idxs_to_remove)

                pressure_segs_clean = pressure_segs_all[patient_id][measurement]['values']
                ecg_clean = pressure_segs_all[patient_id][measurement]['ecg']
                end_idx = np.min([len(x) for x in pressure_segs_clean])
                pressure_segs_same_length = [pressure_segs_clean[II][:end_idx] for II in range(len(pressure_segs_clean))]

                pressure_segs_mean = np.stack(pressure_segs_same_length).mean(axis=0)

                opt_idx_found = False
                unacceptable_idxs = []
                while not opt_idx_found:
                    if len(unacceptable_idxs) > len(pressure_segs_same_length)-1:
                        print(f'none of the ecgs were suitable for {patient_id}')
                        exit()
                    optimal_idx = np.argmin([np.sum((pressure_segs_same_length[II] - 
                        pressure_segs_mean)**2) for II in 
                                         range(len(pressure_segs_same_length)) if 
                                         II not in unacceptable_idxs])

                    # check that the optimal idx has a realistic ecg
                    opt_T = heart_period_all_dict[patient_id][measurement]['T'][optimal_idx]
                    opt_T_vc = heart_period_all_dict[patient_id][measurement]['T_vc'][optimal_idx]
                    opt_T_ac = heart_period_all_dict[patient_id][measurement]['T_ac'][optimal_idx]
                    opt_T_ar = heart_period_all_dict[patient_id][measurement]['T_ar'][optimal_idx]

                    if opt_T_vc < 0.10*opt_T:
                        print('ventricle contraction is unrealistically short',
                              'choosing a new optimal_idx')
                        opt_idx_found = False
                        unacceptable_idxs.append(optimal_idx)
                    elif opt_T_ac < 0.05*opt_T:
                        print('atrial contraction is unrealistically short', 
                              'choosing a new optimal_idx')
                        opt_idx_found = False
                    else:
                        opt_idx_found = True

                pressure_segs_optimal = pressure_segs_clean[optimal_idx]
                ecg_optimal = ecg_clean[optimal_idx]

                fig, ax = plt.subplots()
                [plt.plot(pressure_segs_clean[II], 'k') for II in range(len(pressure_segs_clean))]
                [plt.plot(seg, 'r') for idx, seg in enumerate(pressure_segs) if idx in idxs_to_remove_pressure]
                [plt.plot(seg, 'b') for idx, seg in enumerate(pressure_segs) if idx in idxs_to_remove_ecg]
                plt.plot(pressure_segs_mean, 'g')
                plt.plot(pressure_segs_optimal, 'm')

                plt.savefig(os.path.join(plot_dir, f'pressure_segs_{patient_id}_{measurement}_rem_outliers.eps'))
                plt.close()

                reduced_data_dict[patient_id][measurement] = pressure_segs_optimal
                reduced_data_dict[patient_id][measurement+'_ecg'] = ecg_optimal
                for key in heart_period_all_dict[patient_id][measurement]:
                    reduced_heart_period_dict[patient_id][measurement][key] = \
                        heart_period_all_dict[patient_id][measurement][key][optimal_idx]

                # and plot best ecg
                fig, ax = plt.subplots()
                t = np.linspace(0, reduced_heart_period_dict[patient_id][measurement]['T'], len(ecg_optimal))
                ax.plot(t, ecg_optimal)
                ax.axvline(x=reduced_heart_period_dict[patient_id][measurement]['T_vc'], color='b', label='T_vc_end')
                ax.axvline(x=(reduced_heart_period_dict[patient_id][measurement]['T_vc'] +
                              reduced_heart_period_dict[patient_id][measurement]['T_vr']),
                           color='r', label='T_vr_end')
                ax.axvline(x=reduced_heart_period_dict[patient_id][measurement]['t_astart'], color='m', label='t_astart')
                ax.axvline(x=(reduced_heart_period_dict[patient_id][measurement]['t_astart'] +
                              reduced_heart_period_dict[patient_id][measurement]['T_ac'] +
                              reduced_heart_period_dict[patient_id][measurement]['T_ar']), color='g', label='T_ar_end')
                ax.legend()
                plt.savefig(os.path.join(plot_dir, f'ecg_check_{patient_id}_{measurement}.pdf'))
                plt.close()

            elif measurement in nice_names_tomtec:
                reduced_data_dict[patient_id][measurement] = data_dict[patient_id][measurement][0]
            else:
                continue # don't add extra measurements like 'II'

    # get the mean heart rate and scale all of the traces to that heart rate
    # TODO I should be making sure that the heart_periods are well matched between optimal traces
    heart_period_means = {}
    new_sample_rate = 100

    for patient_id in reduced_heart_period_dict.keys():
        heart_period_means[patient_id] = {}
        fig, ax = plt.subplots()
        period_sums = {}
        for measurement in reduced_heart_period_dict[patient_id].keys():
            for T_type in reduced_heart_period_dict[patient_id][measurement].keys():
                if T_type not in period_sums.keys():
                    period_sums[T_type] = 0
                period = reduced_heart_period_dict[patient_id][measurement][T_type]
                period_sums[T_type] += period

        for T_type in reduced_heart_period_dict[patient_id][measurement].keys():
            heart_period_means[patient_id][T_type] = period_sums[T_type] / \
                                                     len(reduced_heart_period_dict[patient_id])

        n_new_values_old = 9999
        for measurement in reduced_heart_period_dict[patient_id].keys():
            # scale data out to mean heart_period and resample it at 100 Hz
            n_values = len(reduced_data_dict[patient_id][measurement])
            max_time = n_values/sample_rate
            # this switches to 100 Hz then scales to new heart rate
            n_new_values = int(np.floor(max_time * new_sample_rate*reduced_heart_period_dict[patient_id][measurement]['T'] /
                                        heart_period_means[patient_id]['T']))
            # calculate the min length of all series so we can cut off the end to make them the same length
            min_length = min(n_new_values, n_new_values_old)
            n_new_values_old = n_new_values
            new_idxs = [II*sample_rate/new_sample_rate for II in range(n_new_values)]
            reduced_data_dict[patient_id][measurement] = np.interp(new_idxs, [II for II in range(n_values)],
                                                                   reduced_data_dict[patient_id][measurement])
            plt.plot(reduced_data_dict[patient_id][measurement])

        # make all series the same length
        for measurement in reduced_heart_period_dict[patient_id].keys():
            reduced_data_dict[patient_id][measurement] = reduced_data_dict[patient_id][measurement][:min_length]

        plt.savefig(os.path.join(plot_dir, f'pressure_check_{patient_id}_{measurement}.eps'))
        plt.close()

        for key in heart_period_means[patient_id].keys():
            reduced_data_dict[patient_id][key] = heart_period_means[patient_id][key]

    print(reduced_data_dict)
    reduced_df = pd.DataFrame(reduced_data_dict)
    # save invasive data to pickle
    reduced_df.to_pickle(os.path.join(save_path_combined, 'reduced_data.pkl'))

    mmHg_to_Pa = 133.332
    ml_per_s_to_m3_per_s = 1e-6
    mm_per_s_to_m_per_s = 1e-3
    ml_to_m3 = 1e-6
    mm2_to_m2 = 1e-6
    cm2_to_m2 = 1e-4
    mm_to_m = 1e-3
    no_conv = 1.0
    mm_to_m_diam_to_radius = 1e-3/2

    constants_of_interest = [('US.CA.LVOT.DIAM:BMODE', 'r_lvot',
                            'metre', mm_to_m_diam_to_radius, 'biobeat_measurement'),
                              ('US.CA.AO.DIAM_STJ:BMODE', 'r_ao',
                            'metre', mm_to_m_diam_to_radius, 'biobeat_measurement'),
                              ('US.CA.AV.AVA:CONT_VTI:ANTFLOW:DOPPLER', 'A_aov',
                            'm2', cm2_to_m2, 'biobeat_measurement'),
                              ('T', 'T',
                            'second', no_conv, 'biobeat_measurement'),
                             ('T_vc', 'T_vc',
                            'second', no_conv, 'biobeat_measurement'),
                             ('T_vr', 'T_vr',
                              'second', no_conv, 'biobeat_measurement'),
                             ('T_ac', 'T_ac',
                              'second', no_conv, 'biobeat_measurement'),
                             ('T_ar', 'T_ar',
                              'second', no_conv, 'biobeat_measurement'),
                             ('t_astart', 't_astart',
                            'second', no_conv, 'biobeat_measurement')]

    # Now write to json for the
    variables_of_interest = [('US.CA.LA.VOLES:BP:BMODE', 'heart/q_la', 'max', 
                                   'm3', ml_to_m3, 1.0, 21.3, 'q_{la}'),# 1.0),
                                  ('US.CA.LA.VOLED:CPA2D', 'heart/q_la', 'min', 
                                    'm3', ml_to_m3, 1.0, 23.5, 'q_{la}'),#2.0),
                                  ('US.CA.LV.ESV:LVA4D', 'heart/q_lv', 'min', 
                                      'm3', ml_to_m3, 1.0, 18.4, 'q_{lv}'),#2.0),
                                  ('US.CA.LV.EDV:LVA4D', 'heart/q_lv', 'max', 
                                      'm3', ml_to_m3, 1.0, 12.6, 'q_{lv}'),#3.0),
                                  ('US.CA.RA.VOLES:SP:BMODE:A4C', 'heart/q_ra', 
                                      'max', 'm3', ml_to_m3, 1.0, 10.0, 'q_{ra}'),#
                                  ('US.CA.LVOT.VMAX:ANTFLOW:DOPPLER', 'heart/vel_aov', 
                                   'max', 'm_per_s', no_conv, 1.0, 50.0, 'vel_{aov}'), # TODO get the uncertainty
                                  ('LV', 'heart/u_lv', 
                                   ['series', 'max'], 'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'P_{lv}'), # TODO get the uncertainty
                                  ('AO', 'aortic_root/u', 
                                   ['series', 'min', 'max'], 'J_per_m3', mmHg_to_Pa, 1.0, 5.0, 'P_{ao}')] # TODO get the uncertainty
                                  # ('US.CA.LVOT.VMAX:ANTFLOW:DOPPLER', 'heart/vel_aov', 'max', 'alg', 'm_per_s', no_conv, 1.0),#1.0),
                                  # ('US.CA.MV.VMAX_E:ANTFLOW:DOPPLER', 'heart/vel_miv', 'max', 'alg','m_per_s', no_conv, 1.0)]#1.0)]
                                  # ('US.CA.MV.VMAX_A:ANTFLOW:DOPPLER', 'heart/vel_miv', 'max', 'alg','m_per_s', no_conv, 1.0),
                                  # ('US.CA.TV.VMAX:REGFLOW:DOPPLER', 'heart/vel_trv', 'max', 'alg','m_per_s', no_conv, 1.0)]
                                  # ('US.CA.LVOT.DIAM:BMODE', ),
                                  # ('US.CA.VC.DIAM_INF:BMODE', ),
        # variables of interest for the bp-plus participants
        # bp_variables_of_interest = [('cSys', 'aortic_root/u', 'max', 'alg', 'J_per_m3', mmHg_to_Pa, 3.0),#7.0),
        #                             ('cDia', 'aortic_root/u', 'min', 'alg', 'J_per_m3', mmHg_to_Pa, 2.0),#4.0),
        #                             ('cMap', 'aortic_root/u', 'mean', 'alg', 'J_per_m3', mmHg_to_Pa, 1.0)]#4.0)]
        # variables_of_interest = [('Sys', 'brachial_L82/u', 'max', 'alg', 'J_per_m3', mmHg_to_Pa, 1.0, 1.5, 'p_{BR_{LT}}'),#7.0),
        #                             ('Dia', 'brachial_L82/u', 'min', 'alg', 'J_per_m3', mmHg_to_Pa, 1.0, 4.0, 'p_{BR_{LT}}'),#4.0),
        #                             ('Map', 'brachial_L82/u', 'mean', 'alg', 'J_per_m3', mmHg_to_Pa, 1.0, 2.25, 'p_{BR_{LT}}')]#4.0)]
    write_to_json_file(reduced_data_dict, variables_of_interest, constants_of_interest,
                       save_path_combined, sample_rate=100.0)

