import os
import json

# Set up directories
cwd = os.getcwd()
file_dir =  os.path.dirname(cwd)
shadow_data_dir = os.path.join(file_dir, 'shadow_data')
random_pauli_shadow_data_dir = os.path.join(shadow_data_dir, 'random_pauli')
random_clifford_shadow_data_dir = os.path.join(shadow_data_dir, 'random_clifford')
big_random_pauli_shadow_data_dir = os.path.join(shadow_data_dir, 'pauli_10000')
big_random_clifford_shadow_data_dir = os.path.join(shadow_data_dir, 'clifford_10000')
vqe_data_dir = os.path.join(cwd, 'vqe_data')
vqe_1000_data_dir = os.path.join(vqe_data_dir, 'shots_1000')
vqe_10000_data_dir = os.path.join(vqe_data_dir, 'shots_10000')

clifford_shadow_list = os.listdir(random_clifford_shadow_data_dir)
pauli_shadow_list = os.listdir(random_pauli_shadow_data_dir)
big_clifford_shadow_list = os.listdir(big_random_clifford_shadow_data_dir)
big_pauli_shadow_list = os.listdir(big_random_pauli_shadow_data_dir)

results = {"clifford_shadows" : {}, "pauli_shadows" : {}, "big_clifford_shadows" : {}, "big_pauli_shadows" : {}, "vqe_1000" : {}, "vqe_10000" : {}}

print("Clifford Shadows")
for cliff_file in clifford_shadow_list:
    print(cliff_file)
    with open(os.path.join(random_clifford_shadow_data_dir, cliff_file), 'rb') as infile:
        data_dict = json.load(infile)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"], "gs_nrg_tap_cst_up" : data_dict["gs_nrg_tap_cst_up"]}
    results["clifford_shadows"][cliff_file] = important_data


print("Pauli Shadows")
for pauli_file in pauli_shadow_list:
    print(pauli_file)
    with open(os.path.join(random_pauli_shadow_data_dir, pauli_file), 'rb') as infile:
        data_dict = json.load(infile)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"], "gs_nrg_tap_cst_up" : data_dict["gs_nrg_tap_cst_up"]}
    results["pauli_shadows"][pauli_file] = important_data

print("Big Clifford Shadows")
for cliff_file in big_clifford_shadow_list:
    print(cliff_file)
    with open(os.path.join(big_random_clifford_shadow_data_dir, cliff_file), 'rb') as infile:
        data_dict = json.load(infile)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"], "gs_nrg_tap_cst_up" : data_dict["gs_nrg_tap_cst_up"]}
    results["big_clifford_shadows"][cliff_file] = important_data

print("Big Pauli Shadows")
for pauli_file in big_pauli_shadow_list:
    print(pauli_file)
    with open(os.path.join(big_random_pauli_shadow_data_dir, pauli_file), 'rb') as infile:
        data_dict = json.load(infile)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"], "gs_nrg_tap_cst_up" : data_dict["gs_nrg_tap_cst_up"]}
    results["big_pauli_shadows"][pauli_file] = important_data

print("VQE 1000")
for vqe_file in vqe_1000_data_dir:
    print(vqe_file)
    with open(os.path.join(vqe_1000_data_dir, vqe_file), 'rb') as infile:
        data_dict = json.load(infile)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_basic_vqe" : data_dict["gs_nrg_basic_vqe"]}
    results["vqe_1000"][vqe_file] = important_data

print("VQE 10000")
for vqe_file in vqe_10000_data_dir:
    print(vqe_file)
    with open(os.path.join(vqe_10000_data_dir, vqe_file), 'rb') as infile:
        data_dict = json.load(infile)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_basic_vqe" : data_dict["gs_nrg_basic_vqe"]}
    results["vqe_10000"][vqe_file] = important_data


with open(os.path.join(cwd, 'results.json'), 'w') as file:
        json.dump(results, file)