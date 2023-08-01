import os
import json
import numpy as np

# Set up directories
cwd = os.getcwd()
file_dir =  os.path.dirname(cwd)
shadow_data_dir = os.path.join(file_dir, 'shadow_data')
random_pauli_shadow_data_dir = os.path.join(shadow_data_dir, 'random_pauli')
random_clifford_shadow_data_dir = os.path.join(shadow_data_dir, 'random_clifford')
big_random_pauli_shadow_data_dir = os.path.join(shadow_data_dir, 'pauli_10000')
big_random_clifford_shadow_data_dir = os.path.join(shadow_data_dir, 'clifford_10000')
vqe_data_dir = os.path.join(cwd, 'vqe_data')
basic_vqe_data_dir = os.path.join(vqe_data_dir, 'basic')
unitary_partitioning_vqe_data_dir = os.path.join(vqe_data_dir, 'unitary_partitioning')
exact_data_dir = os.path.join(cwd, 'exact_data')

clifford_shadow_list = os.listdir(random_clifford_shadow_data_dir)
pauli_shadow_list = os.listdir(random_pauli_shadow_data_dir)
big_clifford_shadow_list = os.listdir(big_random_clifford_shadow_data_dir)
big_pauli_shadow_list = os.listdir(big_random_pauli_shadow_data_dir)

# results = {"clifford_shadows" : {}, "pauli_shadows" : {}, "big_clifford_shadows" : {}, "big_pauli_shadows" : {}}

# print("Clifford Shadows")
# for cliff_file in clifford_shadow_list:
#     print(cliff_file)
#     with open(os.path.join(random_clifford_shadow_data_dir, cliff_file), 'rb') as infile:
#         data_dict = json.load(infile)
#     important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"], "gs_nrg_tap_cst_up" : data_dict["gs_nrg_tap_cst_up"]}
#     # print("gs nrg = ", data_dict["gs_nrg_tap"])
#     # print("gs nrg cst = ", data_dict["gs_nrg_tap_cst"])
#     # print("gs nrg cst up = ", data_dict["gs_nrg_tap_cst_up"])
#     results["clifford_shadows"][cliff_file] = important_data


# print("Pauli Shadows")
# for pauli_file in pauli_shadow_list:
#     print(pauli_file)
#     with open(os.path.join(random_pauli_shadow_data_dir, pauli_file), 'rb') as infile:
#         data_dict = json.load(infile)
#     important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"], "gs_nrg_tap_cst_up" : data_dict["gs_nrg_tap_cst_up"]}
#     # print("gs nrg = ", data_dict["gs_nrg_tap"])
#     # print("gs nrg cst = ", data_dict["gs_nrg_tap_cst"])
#     # print("gs nrg cst up = ", data_dict["gs_nrg_tap_cst_up"])
#     results["pauli_shadows"][pauli_file] = important_data

# print("Big Clifford Shadows")
# for cliff_file in big_clifford_shadow_list:
#     print(cliff_file)
#     with open(os.path.join(big_random_clifford_shadow_data_dir, cliff_file), 'rb') as infile:
#         data_dict = json.load(infile)
#     important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"], "gs_nrg_tap_cst_up" : data_dict["gs_nrg_tap_cst_up"]}
#     # print("gs nrg = ", data_dict["gs_nrg_tap"])
#     # print("gs nrg cst = ", data_dict["gs_nrg_tap_cst"])
#     # print("gs nrg cst up = ", data_dict["gs_nrg_tap_cst_up"])
#     results["big_clifford_shadows"][cliff_file] = important_data

# print("Big Pauli Shadows")
# for pauli_file in big_pauli_shadow_list:
#     print(pauli_file)
#     with open(os.path.join(big_random_pauli_shadow_data_dir, pauli_file), 'rb') as infile:
#         data_dict = json.load(infile)
#     important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"], "gs_nrg_tap_cst_up" : data_dict["gs_nrg_tap_cst_up"]}
#     # print("gs nrg = ", data_dict["gs_nrg_tap"])
#     # print("gs nrg cst = ", data_dict["gs_nrg_tap_cst"])
#     # print("gs nrg cst up = ", data_dict["gs_nrg_tap_cst_up"])
#     results["big_pauli_shadows"][pauli_file] = important_data

# with open(os.path.join(cwd, 'results.json'), 'w') as file:
#         json.dump(results, file)

with open(os.path.join(cwd, 'results.json'), 'r') as file:
        data = json.load(file)

clifford_results = data["clifford_shadows"]
pauli_results = data["pauli_shadows"]
big_clifford_results = data["big_clifford_shadows"]
big_pauli_results = data["big_pauli_shadows"]

clifford_cst_rel_error = []
clifford_cst_up_rel_error = []
for molecule in clifford_results.keys():
        exact_energy = complex(clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst_up"])
        clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

clifford_cst_vs_cst_up_rel_error = [clifford_cst_rel_error[i] - clifford_cst_up_rel_error[i] for i in range(len(clifford_cst_up_rel_error))]
# Positive values indicate cst + up is better than cst alone
truth_table = [True if c >= 0 else False for c in clifford_cst_vs_cst_up_rel_error]
print(truth_table)

pauli_cst_rel_error = []
pauli_cst_up_rel_error = []
for molecule in pauli_results.keys():
        exact_energy = complex(pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst_up"])
        pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

pauli_cst_vs_cst_up_rel_error = [pauli_cst_rel_error[i] - pauli_cst_up_rel_error[i] for i in range(len(pauli_cst_up_rel_error))]
# Positive values indicate cst + up is better than cst alone
truth_table = [True if c >= 0 else False for c in pauli_cst_vs_cst_up_rel_error]
print(truth_table)

pauli_cst_vs_clifford_cst_rel_error = [pauli_cst_rel_error[i] - clifford_cst_rel_error[i] for i in range(len(pauli_cst_rel_error))]
# Positive values indicate clifford is better than pauli
truth_table = [True if c >= 0 else False for c in pauli_cst_vs_clifford_cst_rel_error]
print(truth_table)



big_pauli_cst_rel_error = []
big_pauli_cst_up_rel_error = []
for molecule in big_pauli_results.keys():
        exact_energy = complex(big_pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst_up"])
        big_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        big_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

big_pauli_cst_vs_cst_up_rel_error = [big_pauli_cst_rel_error[i] - big_pauli_cst_up_rel_error[i] for i in range(len(big_pauli_cst_up_rel_error))]
# Positive values indicate cst + up is better than cst alone
truth_table = [True if c >= 0 else False for c in big_pauli_cst_vs_cst_up_rel_error]
print(truth_table)

big_clifford_cst_rel_error = []
big_clifford_cst_up_rel_error = []
for molecule in big_clifford_results.keys():
        exact_energy = complex(big_clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst_up"])
        big_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        big_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

big_clifford_cst_vs_cst_up_rel_error = [big_clifford_cst_rel_error[i] - big_clifford_cst_up_rel_error[i] for i in range(len(big_clifford_cst_up_rel_error))]
# Positive values indicate cst + up is better than cst alone
truth_table = [True if c >= 0 else False for c in big_clifford_cst_vs_cst_up_rel_error]
print(truth_table)

big_pauli_cst_vs_big_clifford_cst_rel_error = [big_pauli_cst_rel_error[i] - big_clifford_cst_rel_error[i] for i in range(len(big_pauli_cst_rel_error))]
# Positive values indicate clifford is better than pauli
truth_table = [True if c >= 0 else False for c in big_pauli_cst_vs_big_clifford_cst_rel_error]
print(truth_table)

################################################

jw_molecules = []
bk_molecules = []
for molecule in clifford_results.keys():
        # Get rid of .json
        molecule_str = molecule[:-5]
        # separate JW and BK
        if molecule_str[-1] == "W":
                jw_molecules.append(molecule)
        else:
                bk_molecules.append(molecule)

jw_clifford_cst_rel_error = []
jw_clifford_cst_up_rel_error = []
for molecule in jw_molecules:
        exact_energy = complex(clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst_up"])
        jw_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        jw_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

jw_clifford_cst_vs_cst_up_rel_error = [jw_clifford_cst_rel_error[i] - jw_clifford_cst_up_rel_error[i] for i in range(len(jw_clifford_cst_up_rel_error))]
# Positive values indicate cst + up is better than cst alone
truth_table = [True if c >= 0 else False for c in jw_clifford_cst_vs_cst_up_rel_error]
print(truth_table)

bk_clifford_cst_rel_error = []
bk_clifford_cst_up_rel_error = []
for molecule in bk_molecules:
        exact_energy = complex(clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst_up"])
        bk_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        bk_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

bk_clifford_cst_vs_cst_up_rel_error = [bk_clifford_cst_rel_error[i] - bk_clifford_cst_up_rel_error[i] for i in range(len(bk_clifford_cst_up_rel_error))]
# Positive values indicate cst + up is better than cst alone
truth_table = [True if c >= 0 else False for c in bk_clifford_cst_vs_cst_up_rel_error]
print(truth_table)

jw_pauli_cst_rel_error = []
jw_pauli_cst_up_rel_error = []
for molecule in jw_molecules:
        exact_energy = complex(pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst_up"])
        jw_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        jw_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

jw_pauli_cst_vs_cst_up_rel_error = [jw_pauli_cst_rel_error[i] - jw_pauli_cst_up_rel_error[i] for i in range(len(jw_pauli_cst_up_rel_error))]
# Positive values indicate cst + up is better than cst alone
truth_table = [True if c >= 0 else False for c in jw_pauli_cst_vs_cst_up_rel_error]
print(truth_table)

bk_pauli_cst_rel_error = []
bk_pauli_cst_up_rel_error = []
for molecule in bk_molecules:
        exact_energy = complex(pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst_up"])
        bk_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        bk_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

bk_pauli_cst_vs_cst_up_rel_error = [bk_pauli_cst_rel_error[i] - bk_pauli_cst_up_rel_error[i] for i in range(len(bk_pauli_cst_up_rel_error))]
# Positive values indicate cst + up is better than cst alone
truth_table = [True if c >= 0 else False for c in bk_pauli_cst_vs_cst_up_rel_error]
print(truth_table)

jw_pauli_cst_vs_clifford_cst_rel_error = [jw_pauli_cst_rel_error[i] - jw_clifford_cst_rel_error[i] for i in range(len(jw_pauli_cst_rel_error))]
# Positive values indicate clifford is better than pauli
truth_table = [True if c >= 0 else False for c in jw_pauli_cst_vs_clifford_cst_rel_error]
print(truth_table)

bk_pauli_cst_vs_clifford_cst_rel_error = [bk_pauli_cst_rel_error[i] - bk_clifford_cst_rel_error[i] for i in range(len(bk_pauli_cst_rel_error))]
# Positive values indicate clifford is better than pauli
truth_table = [True if c >= 0 else False for c in bk_pauli_cst_vs_clifford_cst_rel_error]
print(truth_table)

# For small clifford compare JW with BK

# For small pauli compare JW with BK

# For big clifford compare JW with BK

# For big pauli compare JW with BK

# Compare small pauli to small clifford

# Compare big pauli to big clifford

# Compare small pauli to big pauli

# Compare small clifford to big clifford
