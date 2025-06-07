from lungsim_impedance_to_gt_outputs import convert_lungsim_output_to_obs_data_json 
from get_intrabeat_periods import generate_intrabeat_periods_json
from create_gt_Alfred import create_gt_obs_data_Alfred

if __name__ == "__main__":

    patient_nums = [4] # [4, 7]
    # project_dir = '/home/farg967/Documents'
    project_dir = '/hpc/farg967/pulmonary_workspace'

    for patient_num in patient_nums:

        use_CO_from_SV = False # if true calcs cardiac output from stroke volume, 
        # otherwise uses CO from ALL_DATA.xlsx file (TODO make sure this isn't hardcoded)
        convert_lungsim_output_to_obs_data_json(patient_num, "pre", project_dir, use_CO_from_SV)
        convert_lungsim_output_to_obs_data_json(patient_num, "post", project_dir, use_CO_from_SV)

        generate_intrabeat_periods_json(patient_num, project_dir)
        create_gt_obs_data_Alfred(patient_num, project_dir, use_CO_from_SV)
    

