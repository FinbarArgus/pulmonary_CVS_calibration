import pandas as pd
import numpy as np
import time
import json
from matplotlib import pyplot as plt
import os
from get_valve_data import write_to_json_file

def create_gt_obs_data_Alfred(patient_num, data_dir, use_CO_from_SV):

    # TODO the file paths should be inputs to this function
    data_file_path = os.path.join(data_dir, 'pulmonary/Alfred_Echo_Pre_Post_vols_P1_13.csv')
    all_data_file_path = os.path.join(data_dir, f'pulmonary/ALL_DATA.xlsx')
    save_dir_path = os.path.join(data_dir, 'pulmonary/ground_truth_for_CA')
    if not os.path.exists(save_dir_path):
        os.mkdir(save_dir_path)
    
    fs = 180
    
    # read csv file 
    df = pd.read_csv(data_file_path, delimiter=',') 

    df_row = df.iloc[np.where(df['Patient No'] == patient_num)]

    # turn df into dict
    data_dict = {patient_num : {entry: df_row[entry].values[0] for entry in df_row.columns}}

    df = pd.read_excel(all_data_file_path)
    # Set the first column as the index
    df.set_index(df.columns[0], inplace=True)

    
    # Find the indices for the pre-surgery and post-surgery sections
    pre_surgery_index = df.index.get_loc('RIGHT HEART CATHETER DATA - pre-surgery')
    post_surgery_index = df.index.get_loc('RIGHT HEART CATHETER DATA - post-surgery')

    # Extract the section between pre-surgery and post-surgery
    pre_surgery_data = df.iloc[pre_surgery_index:post_surgery_index]

    # Extract the 'PCWP Mean' row
    data_dict[patient_num]["P_ao Sys"] = pre_surgery_data.loc['aortic SP'][patient_num]
    data_dict[patient_num]["P_ao Dia"] = pre_surgery_data.loc['aortic DP'][patient_num]
    data_dict[patient_num]["P_pul Sys"] = pre_surgery_data.loc['PASP'][patient_num]
    data_dict[patient_num]["P_pul Dia"] = pre_surgery_data.loc['PADP'][patient_num]
    data_dict[patient_num]["P_rv Sys"] = pre_surgery_data.loc['RVSP'][patient_num]
    data_dict[patient_num]["P_rv Dia"] = pre_surgery_data.loc['RVDP'][patient_num]
    data_dict[patient_num]["P_pcwp mean"] = pre_surgery_data.loc['PCWP Mean'][patient_num]
    data_dict[patient_num]["CO"] = pre_surgery_data.loc['CO'][patient_num]

    ml_to_m3 = 1e-6
    mmHg_to_Pa = 133.332
    no_conv = 1.0
    l_per_min_to_ml_per_m3 = 1e-3/60

    # TODO get the aortic pressure, rv pressure, and PA pressure from dataset here

    variables_of_interest_pre = [('Pre LV Dia Vol- ml', 'heart/q_lv', 'max', 
                                    'm3', ml_to_m3, 1.0, 10.0, 'q_{lv}'),
                                ('Pre LV Sys Vol- ml', 'heart/q_lv', 'min', 
                                    'm3', ml_to_m3, 1.0, 10.0, 'q_{lv}'),
                                ('P_ao Sys', 'ascending_aorta_B/u', 'max', 
                                    'J_per_m3', mmHg_to_Pa, 3.0, 10.0, 'u_{ao}'),
                                ('P_ao Dia', 'ascending_aorta_B/u', 'min', 
                                    'J_per_m3', mmHg_to_Pa, 3.0, 10.0, 'u_{ao}'),
                                ('P_pul Sys', 'MPA_A/u', 'max', 
                                    'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'u_{pul}'),
                                ('P_pul Dia', 'MPA_A/u', 'min', 
                                    'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'u_{pul}'),
                                ('P_rv Sys', 'heart/u_rv', 'max', 
                                    'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'u_{rv}'),
                                ('P_rv Dia', 'heart/u_rv', 'min', 
                                    'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'u_{rv}'),
                                ('P_pcwp mean', 'heart/u_la', 'mean', 
                                    'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'u_{la}')]

    variables_of_interest_post = [('Post LV Dia Vol- ml', 'heart/q_lv', 'max', 
                                    'm3', ml_to_m3, 1.0, 10.0, 'q_{lv}'),
                                ('Post LV Sys Vol- ml', 'heart/q_lv', 'min', 
                                    'm3', ml_to_m3, 1.0, 10.0, 'q_{lv}'),
                                ('P_ao Sys', 'brachial_L82/u', 'max', 
                                    'J_per_m3', mmHg_to_Pa, 3.0, 10.0, 'u_{ao}'),
                                ('P_ao Dia', 'brachial_L82/u', 'min', 
                                    'J_per_m3', mmHg_to_Pa, 3.0, 10.0, 'u_{ao}'),
                                ('P_pul Sys', 'MPA_A/u', 'max', 
                                    'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'u_{pul}'),
                                ('P_pul Dia', 'MPA_A/u', 'min', 
                                    'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'u_{pul}'),
                                ('P_rv Sys', 'heart/u_rv', 'max', 
                                    'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'u_{rv}'),
                                ('P_rv Dia', 'heart/u_rv', 'min', 
                                    'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'u_{rv}'),
                                ('P_pcwp mean', 'heart/u_la', 'mean', 
                                    'J_per_m3', mmHg_to_Pa, 1.0, 10.0, 'u_{la}')]
                                    
    # add atrial mean pressure for the ROM
    constants_of_interest_pre = []
    constants_of_interest_pre.append(('P_pcwp mean', 'u_out_MPA_A',
                            'J_per_m3', mmHg_to_Pa, 'Alfred_database'))
    # constants_of_interest_pre.append(('P_pcwp mean', 'u_out_LPV_V',
    #                         'J_per_m3', mmHg_to_Pa, 'Alfred_database'))
    # constants_of_interest_pre.append(('P_pcwp mean', 'u_out_RPV_V',
    #                         'J_per_m3', mmHg_to_Pa, 'Alfred_database'))

    
    # get periods info
    periods_file = os.path.join(save_dir_path, f'alfred_periods_constants_{patient_num}.json')
    with open(os.path.join(save_dir_path, periods_file), 'r') as file_2:
        data_2 = json.load(file_2)
    # calculate resistances from patient mean pressure, CO and ADAN flows
    # This requires the period from data_2
    terminal_vessels = ["leg_L", "leg_R", 
                        "arm_L", "arm_R", 
                        "posterior_cerebral_L", "posterior_cerebral_R", 
                        "external_carotid_L", "external_carotid_R", 
                        "middle_cerebral_L","middle_cerebral_R",
                        "anterior_cerebral_L", "anterior_cerebral_R", 
                        "trunk_C"]

    ADAN_flows = [7.291e-6, 7.291e-6, 
                  1.548e-6, 1.548e-6,
                  0.564e-6, 0.564e-6,
                  2.7426e-6, 2.7426e-6,
                  2.456e-6, 2.456e-6,
                  1.55e-6, 1.55e-6,
                  37.28e-6]

    ADAN_total = np.sum(ADAN_flows)
    ADAN_fractions = np.array(ADAN_flows)/ADAN_total
    print(data_2)
    period = [data_2[II]["value"] for II in range(len(data_2)) if 
            data_2[II]["variable_name"] == "T"][0]
    # Estimate mean arterial pressure 
    MAP = data_dict[patient_num]["P_ao Dia"] + 1/3*(data_dict[patient_num]["P_ao Sys"] - data_dict[patient_num]["P_ao Dia"])
    MAP = MAP*mmHg_to_Pa

    SV = data_dict[patient_num]['Pre SV (ml)'] 
    SV = SV*ml_to_m3
    CO = SV/period

    flows = ADAN_fractions*CO
    resistances = MAP/flows

    # and get metabolism for each terminal
    C_O2_a = 0.2
    C_O2_p = np.array([0.155, 0.155,
              0.155, 0.155,  
              0.14, 0.14,  
              0.14, 0.14,  
              0.14, 0.14,  
              0.14, 0.14,  
              0.155])

    M_O2 = -flows*(-C_O2_p + C_O2_a)
    M_CO2 = -M_O2*0.75 # TODO improve this assumption for each tissue type.


    for II in range(len(terminal_vessels)):
        R_T_name = 'R_T_' + terminal_vessels[II] + '_T'
        data_dict[patient_num][R_T_name] = resistances[II]
        constants_of_interest_pre.append((R_T_name, R_T_name,
                            'Js_per_m6', no_conv, 'Alfred_ADAN_calculated'))
        M_O2_name = 'M_O2_' + terminal_vessels[II] + '_GE'
        M_CO2_name = 'M_CO2_' + terminal_vessels[II] + '_GE'
        data_dict[patient_num][M_O2_name] = M_O2[II]
        constants_of_interest_pre.append((M_O2_name, M_O2_name,
                            'm3_per_s', no_conv, 'Alfred_ADAN_Ursino_calculated'))
        data_dict[patient_num][M_CO2_name] = M_CO2[II]
        constants_of_interest_pre.append((M_CO2_name, M_CO2_name,
                            'm3_per_s', no_conv, 'Alfred_ADAN_Ursino_calculated'))

    # add cardiac output here for the ROM
    data_dict[patient_num]["CO (from SV)"] = CO
    print("CO from catheter measurement = ", 
            data_dict[patient_num]["CO"]*l_per_min_to_ml_per_m3, "m3_per_s")
    print("CO from ECHO SV and ECG period = ", 
            data_dict[patient_num]["CO (from SV)"], "m3_per_s")
    if use_CO_from_SV:
        constants_of_interest_pre.append(('CO (from SV)', 'A_0_inlet',
                                'm3_per_s', no_conv, 'Alfred_database'))
        print("using CO from ECHO SV and ECG period for ROM")
    else:
        constants_of_interest_pre.append(('CO', 'A_0_inlet',
                                'm3_per_s', l_per_min_to_ml_per_m3, 'Alfred_database'))
        print("using CO from catheter for ROM")

    # add resistances here
    
    save_name = 'ground_truth_Alfred'
    save_name_pre = save_name + '_pre'
    write_to_json_file(data_dict, variables_of_interest_pre, constants_of_interest_pre,
                    save_dir_path, save_name=save_name_pre, sample_rate=fs)

    # combine periods json constants file with the one just generated for the pre case

    with open(os.path.join(save_dir_path, f'{save_name_pre}_constants_{patient_num}.json'), 'r') as file_1:
        data_1 = json.load(file_1)

    # also add vessel geometry information
    geom_file_pre = os.path.join(save_dir_path, f'ROM_gt/vessel_geom_constants_pre_patient_{patient_num}.json')
    geom_file_post = os.path.join(save_dir_path, f'ROM_gt/vessel_geom_constants_post_patient_{patient_num}.json')
    
    with open(os.path.join(save_dir_path, geom_file_pre), 'r') as file_3_pre:
        data_3_pre = json.load(file_3_pre)
    
    with open(os.path.join(save_dir_path, geom_file_post), 'r') as file_3_post:
        data_3_post = json.load(file_3_post)

    # ths only works because my json data is a list
    combined_data_pre = data_1 + data_2 + data_3_pre
    print(combined_data_pre)
    # post case still has periods from pre, because we assume they are unknown for post
    combined_data_post = data_1 + data_2 + data_3_post

    with open(os.path.join(save_dir_path, f'{save_name}_pre_constants_{patient_num}.json'), 'w') as wf:
        json.dump(combined_data_pre, wf, indent=2)
    
    with open(os.path.join(save_dir_path, f'{save_name}_post_constants_{patient_num}.json'), 'w') as wf:
        json.dump(combined_data_post, wf, indent=2)


if __name__ == '__main__':
        
    if len(sys.argv) == 3:
        patient_num=sys.argv[1]
        data_dir =sys.argv[2]
        
    else:
        print("usage:  python create_gt_Alfred.py patient_num data_dir") 
        exit()

    create_gt_obs_data_Alfred(patient_num, data_dir)
