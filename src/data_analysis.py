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
random_pauli_shadow_data_dir = os.path.join(shadow_data_dir, 'random_pauli')
random_clifford_shadow_data_dir = os.path.join(shadow_data_dir, 'random_clifford')
big_random_pauli_shadow_data_dir = os.path.join(shadow_data_dir, 'pauli_10000')
big_random_clifford_shadow_data_dir = os.path.join(shadow_data_dir, 'clifford_10000')
vqe_data_dir = os.path.join(file_dir, 'vqe_data')
vqe_1000_data_dir = os.path.join(vqe_data_dir, 'shots_1000')
vqe_10000_data_dir = os.path.join(vqe_data_dir, 'shots_10000')

clifford_shadow_list = os.listdir(random_clifford_shadow_data_dir)
pauli_shadow_list = os.listdir(random_pauli_shadow_data_dir)
big_clifford_shadow_list = os.listdir(big_random_clifford_shadow_data_dir)
big_pauli_shadow_list = os.listdir(big_random_pauli_shadow_data_dir)
vqe_1000_list = os.listdir(vqe_1000_data_dir)
vqe_10000_list = os.listdir(vqe_10000_data_dir)

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

print("Big Pauli Shadows")
for pauli_file in big_pauli_shadow_list:
    print(pauli_file)
    with open(os.path.join(big_random_pauli_shadow_data_dir, pauli_file), 'rb') as infile:
        data_dict = json.load(infile)

    with open(os.path.join(ham_data_dir, pauli_file), 'r') as infile:
        ham_file = json.load(infile)
    H = PauliwordOp.from_dictionary(ham_file['hamiltonian'])

    QT = QubitTapering(H)
    hf_state   = QuantumState(np.asarray(ham_file['data']['hf_array'])) # Hartree-Fock state
    hf_energy  = ham_file['data']['calculated_properties']['HF']['energy']
    H_taper   = QT.taper_it(ref_state=hf_state) 
    gs_nrg_tap, gs_psi_tap = exact_gs_energy(H_taper.to_sparse_matrix)
    clique_cover = H_taper.clique_cover(edge_relation='C')

    classical_shadow = ClassicalShadow(QuantumCircuit(gs_psi_tap.n_qubits), [])
    classical_shadow.unitary_ensemble = "pauli"
    classical_shadow.shadows = []
    for s, b in data_dict["shadow"].values():
        classical_shadow.shadows.append((s,b))
    classical_shadow.observables = []
    for clique in clique_cover.values():
        classical_shadow.observables.append(clique.to_sparse_matrix)
    classical_shadow.num_shadows = len(classical_shadow.shadows)

    _, res = classical_shadow.linearPredictions(10)
    gs_nrg_tap_cst_gc = str(sum(res))
    print(data_dict["gs_nrg_tap"], data_dict["gs_nrg_tap_cst"], data_dict["gs_nrg_tap_cst_up"], gs_nrg_tap_cst_gc)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"], "gs_nrg_tap_cst_up" : data_dict["gs_nrg_tap_cst_up"], "gs_nrg_tap_cst_gc" : gs_nrg_tap_cst_gc}

    results["big_pauli_shadows"][pauli_file] = important_data

print("Big Clifford Shadows")
for cliff_file in big_clifford_shadow_list:
    print(cliff_file)
    with open(os.path.join(big_random_clifford_shadow_data_dir, cliff_file), 'rb') as infile:
        data_dict = json.load(infile)

    with open(os.path.join(ham_data_dir, cliff_file), 'r') as infile:
        ham_file = json.load(infile)
    H = PauliwordOp.from_dictionary(ham_file['hamiltonian'])
    QT = QubitTapering(H)
    hf_state   = QuantumState(np.asarray(ham_file['data']['hf_array'])) # Hartree-Fock state
    hf_energy  = ham_file['data']['calculated_properties']['HF']['energy']
    H_taper   = QT.taper_it(ref_state=hf_state) 
    gs_nrg_tap, gs_psi_tap = exact_gs_energy(H_taper.to_sparse_matrix)
    clique_cover = H_taper.clique_cover(edge_relation='C')

    classical_shadow = ClassicalShadow(QuantumCircuit(gs_psi_tap.n_qubits), [])
    classical_shadow.unitary_ensemble = "random clifford"
    classical_shadow.shadows = []
    for s, b in data_dict["shadow"].values():
        classical_shadow.shadows.append((s,b))

    classical_shadow.observables = []
    for clique in clique_cover.values():
        classical_shadow.observables.append(clique.to_sparse_matrix)
    classical_shadow.num_shadows = len(classical_shadow.shadows)
    _, res = classical_shadow.linearPredictions(10)
    gs_nrg_tap_cst_gc = str(sum(res))
    print(data_dict["gs_nrg_tap"], data_dict["gs_nrg_tap_cst"], data_dict["gs_nrg_tap_cst_up"], gs_nrg_tap_cst_gc)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_tap_cst" : data_dict["gs_nrg_tap_cst"], "gs_nrg_tap_cst_up" : data_dict["gs_nrg_tap_cst_up"], "gs_nrg_tap_cst_gc" : gs_nrg_tap_cst_gc}

    results["big_clifford_shadows"][cliff_file] = important_data

print("VQE 1000")
for vqe_file in vqe_1000_list:
    print(vqe_file)
    with open(os.path.join(vqe_1000_data_dir, vqe_file), 'rb') as infile:
        data_dict = json.load(infile)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_basic_vqe" : data_dict["gs_nrg_basic_vqe"]}
    results["vqe_1000"][vqe_file] = important_data

print("VQE 10000")
for vqe_file in vqe_10000_list:
    print(vqe_file)
    with open(os.path.join(vqe_10000_data_dir, vqe_file), 'rb') as infile:
        data_dict = json.load(infile)
    important_data = {"gs_nrg_tap" : data_dict["gs_nrg_tap"], "gs_nrg_basic_vqe" : data_dict["gs_nrg_basic_vqe"]}
    results["vqe_10000"][vqe_file] = important_data


with open(os.path.join(cwd, 'results.json'), 'w') as file:
        json.dump(results, file)