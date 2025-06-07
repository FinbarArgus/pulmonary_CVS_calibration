import os, sys
import pandas as pd
import json
root_dir_path = os.path.join(os.path.dirname(__file__), '../../')
sys.path.append(os.path.join(root_dir_path, 'python_ECG_tools'))
sys.path.append(os.path.join(root_dir_path, 'ultrasound_uncertainty_tools'))
sys.path.append(os.path.join(root_dir_path, 'circulatory_autogen'))

from read_tomtec_worksheets import read_tomtec_worksheets as rt
from haemotools.BB_echo_functions import BBdata
from read_tomtec_worksheets import read_tomtec_worksheets

inv_data_dir = '/mnt/heart-eresearch/Projects/biobeat/data/pressures/invasive'
tomtec_data_dir = '/mnt/heart-eresearch/Projects/biobeat/analyses/echo-reporting/tomtec-worksheets/'

save_path_inv = '/mnt/heart-eresearch/Sandboxes/Finbar/invasive/'
save_path_tomtec = '/mnt/heart-eresearch/Sandboxes/Finbar/tomtec/'

ec = BBdata()
inv_data = ec.data_dict_general(inv_data_dir)
# result = json_df.to_json(orient='records')
# parsed = json.loads(result)
    # json.dump(parsed, wf, indent=2)
    # json.dump(inv_data, wf, indent=2)

inv_df = pd.DataFrame(inv_data)
# save invasive data to pickle
inv_df.to_pickle(os.path.join(save_path_inv, 'inv_data.pkl'))

study_data = rt.read_series_of_tomtec_worksheets_s(tomtec_data_dir, 'BB')

rt.write_data_to_csv_file_s(study_data, os.path.join(save_path_tomtec, 'tomtec_data.csv'))
rt.write_desc_to_csv_file_s(study_data, os.path.join(save_path_tomtec, 'tomtec_desc.csv'))
