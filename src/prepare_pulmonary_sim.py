import os
import sys
import yaml
from change_patient_num import change_patient_num

root_dir = os.path.join(os.path.dirname(__file__), '../../circulatory_autogen') # this assumes circulatory_autogen is in the same directory as lv_pressure_estimation
sys.path.append(os.path.join(root_dir, 'src/scripts'))
user_inputs_dir = os.path.join(root_dir, 'user_run_files')

from read_and_insert_parameters import insert_parameters

def prepare_pulmonary_sim(patient_num, case_type, project_dir):
    change_patient_num(patient_num, case_type, project_dir)

    # TODO change sim times here.

    with open(os.path.join(user_inputs_dir, 'user_inputs.yaml'), 'r') as file:
        inp_data_dict = yaml.load(file, Loader=yaml.FullLoader)

    resources_dir = os.path.join(root_dir, 'resources')
    # overwrite dir paths if set in user_inputs.yaml
    if "resources_dir" in inp_data_dict.keys():
        resources_dir = inp_data_dict['resources_dir']
        
    # TODO should this be an input. 
    if case_type.endswith('pre'):
        case_tshort = 'pre'
    elif case_type.endswith('post'):
        case_tshort = 'post'

    parameters_csv_abs_path = os.path.join(resources_dir, inp_data_dict['input_param_file'])


    constants_path = os.path.join(project_dir, f"data/pulmonary/ground_truth_for_CA/ground_truth_Alfred_{case_tshort}_constants_{patient_num}.json")
    insert_parameters(parameters_csv_abs_path, constants_path)
    if case_type.startswith('coupled'):
        pulmonary_model_path = os.path.join(project_dir, f'physiology_models/pulmonary_CVS_Alfred/patient_{patient_num}/{case_tshort}/' + \
                               f'generated_models/lung_ROM_lung_ROM_lobe_imped_{case_tshort}_patient_{patient_num}_obs_data/lung_ROM_parameters.csv')
        insert_parameters(parameters_csv_abs_path, pulmonary_model_path)

if __name__ == "__main__":
    if len(sys.argv) == 3:
        patient_num=sys.argv[1]
        case_type=sys.argv[2]
        if case_type not in ['pre', 'post', 'coupled_pre', 'coupled_post']:
            print(f'case type must be in [pre, post, coupled_pre, coupled_post], not {case_type}')
            exit()

    else:
        print("usage: python prepare_pulmonary_sim.py patient_num case_type") 
        exit()
    
    # Change this if you want to move the project workspace directory
    project_dir = "/hpc/farg967/pulmonary_workspace"
    print(f'preparing pulmonary simulations in dir {project_dir}')
    prepare_pulmonary_sim(patient_num, case_type, project_dir)


