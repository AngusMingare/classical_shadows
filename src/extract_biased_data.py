import os
import json
from symmer import PauliwordOp, QubitTapering, QuantumState
from tomography.cst import ClassicalShadow
import numpy as np
from qiskit import QuantumCircuit
from symmer.utils import exact_gs_energy

# Set up directories
cwd = os.getcwd()
file_dir =  os.path.dirname(cwd)
symmer_dir = os.path.join(cwd, 'symmer')
test_dir = os.path.join(symmer_dir, 'tests')
ham_data_dir = os.path.join(test_dir, 'hamiltonian_data')
shadow_data_dir = os.path.join(file_dir, 'shadow_data')
big_random_clifford_shadow_data_dir = os.path.join(shadow_data_dir, 'clifford_10000')
big_biased_clifford_shadow_data_dir = os.path.join(shadow_data_dir, 'biased_clifford')

temp_big_clifford_shadow_list = os.listdir(big_random_clifford_shadow_data_dir)
big_biased_clifford_shadow_list = os.listdir(big_biased_clifford_shadow_data_dir)
big_biased_clifford_shadow_list.remove("H2_6-31G_SINGLET_BK.json")
big_biased_clifford_shadow_list.remove("H2_3-21G_SINGLET_BK.json")
big_biased_clifford_shadow_list.remove("H2_6-31G_SINGLET_JW.json")
big_biased_clifford_shadow_list.remove("H2_3-21G_SINGLET_JW.json")
big_biased_clifford_shadow_list.remove("HeH+_3-21G_SINGLET_BK.json")
big_biased_clifford_shadow_list.remove("HeH+_3-21G_SINGLET_JW.json")

big_clifford_shadow_list = [x for x in temp_big_clifford_shadow_list if x in big_biased_clifford_shadow_list]

results = {"big_clifford_shadows" : {}, "big_biased_clifford_shadows" : {}}

print("Big Clifford Shadows")
for cliff_file in big_clifford_shadow_list:
    print(cliff_file)
    with open(os.path.join(big_random_clifford_shadow_data_dir, cliff_file), 'rb') as infile:
        data_dict = json.load(infile)

    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"]}

    results["big_clifford_shadows"][cliff_file] = important_data

print("Big Biased Clifford Shadows")
for cliff_file in big_biased_clifford_shadow_list:
    print(cliff_file)
    with open(os.path.join(big_biased_clifford_shadow_data_dir, cliff_file), 'rb') as infile:
        data_dict = json.load(infile)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"]}

    results["big_biased_clifford_shadows"][cliff_file] = important_data

with open(os.path.join(cwd, 'biased_results.json'), 'w') as file:
        json.dump(results, file)