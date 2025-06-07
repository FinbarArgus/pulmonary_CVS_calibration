from lungsim_impedance_to_gt_outputs import convert_lungsim_output_to_obs_data_json 
from get_intrabeat_periods import generate_intrabeat_periods_json
from create_gt_Alfred import create_gt_obs_data_Alfred
import shutil
import os

def copy_patient_resources(patient_num):
    srcs = ["../patient_models/patient_generic/pre/resources",
            "../patient_models/patient_generic/post/resources",
            "../patient_models/patient_generic/coupled_pre/resources"]
            # "../patient_models/patient_generic/coupled_post/resources"] # todo setup code for coupling calibrated post ROM to CVS
    dests = [f"../patient_models/patient_{patient_num}/pre/resources",
            f"../patient_models/patient_{patient_num}/post/resources",
            f"../patient_models/patient_{patient_num}/coupled_pre/resources"]
            # f"../patient_models/patient_{patient_num}/coupled_post/resources"] # todo setup code for coupling calibrated post ROM to CVS
    
    for src, dest in zip(srcs, dests):
        # Make sure the destination parent directory exists
        os.makedirs(os.path.dirname(os.path.dirname(dest)), exist_ok=True)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        
        # Copy directory (dir must not already exist)
        shutil.copytree(src, dest, dirs_exist_ok=True)
    

if __name__ == "__main__":

    patient_nums = [4] # [4, 7]
    data_dir = '/eresearch/abi-12l-ep1/farg967/farg967_sandbox/data'

    for patient_num in patient_nums:
        # copy patient_generic/resources to patient_num/resources
        copy_patient_resources(patient_num)

        use_CO_from_SV = False # if true calcs cardiac output from stroke volume, 
        # otherwise uses CO from ALL_DATA.xlsx file (TODO make sure this isn't hardcoded)
        convert_lungsim_output_to_obs_data_json(patient_num, "pre", data_dir, use_CO_from_SV)
        convert_lungsim_output_to_obs_data_json(patient_num, "post", data_dir, use_CO_from_SV)

        generate_intrabeat_periods_json(patient_num, data_dir)
        create_gt_obs_data_Alfred(patient_num, data_dir, use_CO_from_SV)
    

