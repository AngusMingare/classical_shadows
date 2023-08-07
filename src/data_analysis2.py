import numpy as np
import json
import os 

# Set up directories
cwd = os.getcwd()
file_dir =  os.path.dirname(cwd)

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

