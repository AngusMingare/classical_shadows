import os
import json
import numpy as np
from symmer import PauliwordOp, QuantumState, QubitTapering
from symmer.operators import AntiCommutingOp
from symmer.utils import exact_gs_energy
from symmer.evolution.exponentiation import *
from symmer.evolution.decomposition import *
from tomography.cst import ClassicalShadow
from qiskit import QuantumCircuit
from scipy.sparse import csr_array

# Set up directories
cwd = os.getcwd()
file_dir =  os.path.dirname(cwd)
symmer_dir = os.path.join(cwd, 'src/symmer')
test_dir = os.path.join(symmer_dir, 'tests')
ham_data_dir = os.path.join(test_dir, 'hamiltonian_data')
shadow_data_dir = os.path.join(cwd, 'shadow_data')
biased_clifford_data_dir = os.path.join(shadow_data_dir, 'biased_clifford')
globally_biased_pauli_data_dir = os.path.join(shadow_data_dir, 'globally_biased_pauli')

# Hamiltonian Data
ham_list = os.listdir(ham_data_dir)
num_qubits = []
for ham in ham_list:
    with open(os.path.join(ham_data_dir, ham), 'r') as infile:
        data_dict = json.load(infile)
    H = PauliwordOp.from_dictionary(data_dict['hamiltonian'])
    num_qubits.append(H.n_qubits)
sorted_list = sorted(zip(ham_list, num_qubits), key = lambda x: x[1])
ham_list = [ham for ham,_ in sorted_list]

total_measurement_budget = 10000

for filename in ham_list[:20]:

    print(filename)

    with open(os.path.join(ham_data_dir, filename), 'r') as infile:
        data_dict = json.load(infile)
    H = PauliwordOp.from_dictionary(data_dict['hamiltonian'])

    # Tapering
    QT = QubitTapering(H)
    hf_state   = QuantumState(np.asarray(data_dict['data']['hf_array'])) # Hartree-Fock state
    hf_energy  = data_dict['data']['calculated_properties']['HF']['energy']
    H_taper   = QT.taper_it(ref_state=hf_state) 

    # Exact Energy
    gs_nrg_tap, gs_psi_tap = exact_gs_energy(H_taper.to_sparse_matrix)
    print("gs_nrg_tap = ", gs_nrg_tap)

    # Get qiskit circuit preparing gs_psi_tap
    sp_circ = QuantumCircuit(gs_psi_tap.n_qubits)
    gs_psi_tap_dict = gs_psi_tap.to_dictionary
    gs_psi_tap_array = np.array([0]*(2**gs_psi_tap.n_qubits), dtype=complex)
    def bin_to_int(b):
        l = len(b)
        val = 0
        for a in range(l):
            val += int(b[a]) * (2**(l-1-a))
        return val
    for p, w in zip(list(gs_psi_tap_dict.keys()), list(gs_psi_tap_dict.values())):
        val = bin_to_int(p)
        gs_psi_tap_array[val] = w 
    sp_circ.prepare_state(gs_psi_tap_array)
    gs_psi_tap_array = csr_array(gs_psi_tap_array).transpose()

    clique_cover = H_taper.clique_cover(edge_relation='AC')
    commuting_clique_cover = H_taper.clique_cover(edge_relation='C')

  
    # # Estimate gs energy given gs_psi_tap and H_taper using classical shadows
    sp_circ_copy = sp_circ.copy()
    H_taper_dict = H_taper.to_dictionary
    obs = list(H_taper_dict.keys())
    classical_shadow = ClassicalShadow(sp_circ_copy, obs, H_taper)
    classical_shadow.createClassicalShadows(unitary_ensemble="biased clifford", num_shadows=total_measurement_budget)

    classical_shadow.observables = [PauliwordOp.from_dictionary({o : w}).to_sparse_matrix for o,w in zip(list(H_taper_dict.keys()), list(H_taper_dict.values()))]
    obs, results = classical_shadow.linearPredictions(10)
    gs_nrg_tap_cst = 0
    for w, exp in zip(list(H_taper_dict.values()), results):
        gs_nrg_tap_cst += exp 
    print("gs_nrg_tap_cst = ", str(gs_nrg_tap_cst))

    # classical_shadow.total_measurement_budget = 1000
    # estimate = classical_shadow.OGMPauliEstimator("importance sampling")
    # print("estimate = ", str(estimate))

    # Estimate gs energy given gs_psi_tap and H_taper using unitary partitioning and classical shadows
    classical_shadow.observables = []
    weights = []
    for clique in clique_cover.values():
        ac_op = AntiCommutingOp.from_PauliwordOp(clique)
        pop,rots,w,op = ac_op.unitary_partitioning()
        classical_shadow.observables.append(op.to_sparse_matrix)
        weights.append(w)
    _, results = classical_shadow.linearPredictions(10)
    gs_nrg_tap_cst_up = 0
    for w, exp in zip(weights, results):
        gs_nrg_tap_cst_up += w * exp 
    print("gs_nrg_tap_cst_up = ", gs_nrg_tap_cst_up)

    # Estimate gs energy given gs_psi_tap and H_taper using commuting grouping and classical shadows
    classical_shadow.observables = []
    for clique in commuting_clique_cover.values():
        classical_shadow.observables.append(clique.to_sparse_matrix)
    _, res = classical_shadow.linearPredictions(10)
    gs_nrg_tap_cst_gc = str(sum(res))
    print("gs_nrg_tap_cst_gc = ", gs_nrg_tap_cst_gc)

    with open(os.path.join(biased_clifford_data_dir, filename), 'w') as file:
        shadow_data = {}
        for i in range(len(classical_shadow.shadows)):
            matrix, bitstring = classical_shadow.shadows[i]
            matrix = matrix.qasm()
            shadow_data[i] = (matrix, bitstring)
        data = {"num_shadows" : total_measurement_budget,
                "gs_nrg_tap" : str(gs_nrg_tap),
                "gs_nrg_tap_cst" : str(gs_nrg_tap_cst),
                "gs_nrg_tap_cst_up" : str(gs_nrg_tap_cst_up),
                "gs_nrg_tap_cst_gc" : str(gs_nrg_tap_cst_gc),
                "shadow" : shadow_data
                }
        json.dump(data, file)