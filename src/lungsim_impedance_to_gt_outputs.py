import json
import pandas as pd
import sys
import os

def convert_lungsim_output_to_obs_data_json(patient_num, pre_or_post, data_dir, use_CO_from_SV):

    data_file_path = os.path.join(data_dir, f'pulmonary/lobe_impedances/{patient_num}/{pre_or_post}/lobe_imped.json')
    vessel_array_path = os.path.join(os.path.dirname(__file__), f"../patient_models/patient_{patient_num}/{pre_or_post}/resources/lung_ROM_vessel_array.csv")
    save_file_path = os.path.join(data_dir, f'pulmonary/ground_truth_for_CA/ROM_gt/lung_ROM_lobe_imped_{pre_or_post}_patient_{patient_num}_obs_data.json')
    constants_save_file_path = os.path.join(data_dir, f'pulmonary/ground_truth_for_CA/ROM_gt/vessel_geom_constants_{pre_or_post}_patient_{patient_num}.json')
    all_data_file_path = os.path.join(data_dir, f'pulmonary/ALL_DATA.xlsx')

    with open(data_file_path, 'r') as file:
        data = json.load(file)

    vessel_array = pd.read_csv(vessel_array_path, index_col=0)
    # dyne_s_per_cm5_to_J_s_per_m6 
    conversion = 1e5
    l_per_min_to_ml_per_m3 = 1e-3/60
    mmHg_to_Pa = 133.32

    # P_pcwp_mean = 12*mmHg_to_Pa # TODO get this from the ALL_DATA file
    
    df = pd.read_excel(all_data_file_path)
    # Set the first column as the index
    df.set_index(df.columns[0], inplace=True)

    
    # Find the indices for the pre-surgery and post-surgery sections
    pre_surgery_index = df.index.get_loc('RIGHT HEART CATHETER DATA - pre-surgery')
    post_surgery_index = df.index.get_loc('RIGHT HEART CATHETER DATA - post-surgery')

    # Extract the section between pre-surgery and post-surgery
    pre_surgery_data = df.iloc[pre_surgery_index:post_surgery_index]

    # Extract the 'PCWP Mean' row
    P_pcwp_mean = pre_surgery_data.loc['PCWP Mean'][patient_num]*mmHg_to_Pa
    MPA_mean_pressure = pre_surgery_data.loc['PA Mean'][patient_num]*mmHg_to_Pa 

    mu = 0.004

    full_dict = {}
    entry_list = []
    constant_list = []
    terminal_names = []
    main_arteries = ["MPA_A"]
    for II in range(len(data["vessel_names"])):

        entry = {}
        entry["variable"] = data["vessel_names"][II]

        if entry["variable"] not in main_arteries:
            continue
        
        resistance_entry = {}
        resistance_entry["variable_name"] = f'R_{entry["variable"]}'
        resistance_entry["units"] = 'Js_per_m6'
        resistance_entry["data_reference"] = 'Alfred_database'
        resistance_entry["value"] = conversion*data["impedance"][data["vessel_names"][II]][0]
        constant_list.append(resistance_entry)
        print("resistance", resistance_entry["value"])

        entry["data_type"] = "frequency"
        entry["operation"] = "division"
        input_vessel = vessel_array["inp_vessels"][data["vessel_names"][II]].strip()
        print(input_vessel)
        BC_type = vessel_array["BC_type"][data["vessel_names"][II]].strip()
        entry["operands"] = [f'{data["vessel_names"][II]}/u',
                             f'{input_vessel}/v']

        updated_imp = [val*conversion for idx, val in enumerate(data["impedance"][data["vessel_names"][II]]) \
                if not 0 < idx < 10]
        updated_std = [val/10 for val in updated_imp]
        updated_phase = [val for idx, val in enumerate(data["phase"][data["vessel_names"][II]]) \
                if not 0 < idx < 10]
        updated_freq = [val for idx, val in enumerate(data["frequency"]) \
                if not 0 < idx < 10]
        entry["unit"] = "Js/m^6" # data["impedance"]["unit"]
        entry["obs_type"] = "frequency"
        entry["value"] = updated_imp
        entry["std"] = updated_std
        entry["frequencies"] = updated_freq
        entry["phase"] = updated_phase

        entry["weight"] = [0.0 for val in updated_imp]
        entry["phase_weight"] = [0.0 for val in updated_imp]
        for idx in range(len(updated_imp)):
            if idx == 0 :
                # entry["weight"][idx] = 0.0# big because total resistance is most important
                entry["weight"][0] = 0.0 # zero because the total resistance isn't free
            elif idx < 10:
                entry["weight"][idx] = 0
            elif idx <30:
                entry["weight"][idx] = 10
            elif idx <40:
                entry["weight"][idx] = 0.1
            elif idx <50:
                entry["weight"][idx] = 0.01
            elif idx <60:
                entry["weight"][idx] = 0.01
            else:
                entry["weight"][idx] = 0.001


            if idx<10:
                entry["phase_weight"][idx] = 0.0 # entries less than HR period will always be zero
            elif idx<20:
                entry["phase_weight"][idx] = 0.5 
            elif idx<30:
                entry["phase_weight"][idx] = 0.1
            else:
                entry["phase_weight"][idx] = 0.0

        # get the mean flow for this vessel
        if data['mean flow']['unit'] in ['mm3/s', 'mm^3/s']:
            flow_conversion = 1e-9
        else:
            print(f'flow unit of {data["mean flow"]["unit"]} is unknown')
            exit()
        mean_flow = data['mean flow'][data["vessel_names"][II]][0]*flow_conversion

        # add on the downstream resistance past the left atrium ( a ficticious resistance)
        # this is needed becuase in CA we just do the freq domain conversion on Pressure/flow
        entry["value"][0] = entry["value"][0] + P_pcwp_mean/mean_flow

        entry_list.append(entry)

        # create entry for fitting the mean flow
        flow_entry = {}
        flow_entry["variable"] = data["vessel_names"][II] + '/v'
        flow_entry["data_type"] = "constant"
        flow_entry["unit"] = "m3/s" # data["impedance"]["unit"]
        flow_entry["obs_type"] = "mean"
        flow_entry["value"] = mean_flow
        flow_entry["std"] = 0.1*mean_flow
        flow_entry["weight"] = 0.1
        entry_list.append(flow_entry)

        # if post get MPA input resistance to calc mean pressure and overwrite 
        # measured mean pressure
        if pre_or_post == 'post':
            if entry["variable"] == "MPA_A":
                MPA_resistance = entry["value"][0]
                MPA_mean_pressure = MPA_resistance*mean_flow # TODO should this be 
                                                             # The MPA mean pressure 
                                                             # from ALL_DATA

    # add entries for MPA pressure
    # I think this makes sure the model converges quickly. Not super sure
    entry = {}
    entry["variable"] = "MPA_A/u"
    entry["data_type"] = "constant"
    entry["unit"] = "J/m3" # data["impedance"]["unit"]
    entry["obs_type"] = "mean"
    entry["value"] = MPA_mean_pressure
    entry["std"] = 0.1*MPA_mean_pressure
    entry["weight"] = 3
    entry_list.append(entry)

    full_dict["data_item"] = entry_list 

    with open(save_file_path, 'w') as wf: 
        json.dump(full_dict, wf, indent=2)

    with open(constants_save_file_path, 'w') as wf:
        json.dump(constant_list, wf)

if __name__ == "__main__":

    if len(sys.argv) == 4:
        patient_num=sys.argv[1]
        pre_or_post = sys.argv[2]
        data_dir = sys.argv[3]
        if pre_or_post not in ['pre', 'post']:
            print(f'pre_or_post must be "pre" or "post", not {pre_or_post}')
            exit()
        
    else:
        print("usage:  python lungsim_impedance_to_gt_output.py patient_num pre_or_post data_dir") 
        exit()

    convert_lungsim_output_to_obs_data_json(patient_num, pre_or_post, data_dir)



