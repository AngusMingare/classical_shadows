import os
import json
import numpy as np
from symmer import PauliwordOp, QuantumState, QubitTapering
from symmer.operators import AntiCommutingOp
from symmer.utils import exact_gs_energy
from symmer.evolution.exponentiation import *
from symmer.evolution.decomposition import *
from vqe.measurements import measureObservable
from tomography.cst import ClassicalShadow
from qiskit import QuantumCircuit
from scipy.sparse import csr_array, csr_matrix, kron
from functools import reduce

# Set up directories
cwd = os.getcwd()
file_dir =  os.path.dirname(cwd)
symmer_dir = os.path.join(cwd, 'src/symmer')
test_dir = os.path.join(symmer_dir, 'tests')
ham_data_dir = os.path.join(test_dir, 'hamiltonian_data')
shadow_data_dir = os.path.join(cwd, 'shadow_data')
random_pauli_shadow_data_dir = os.path.join(shadow_data_dir, 'random_pauli')
random_clifford_shadow_data_dir = os.path.join(shadow_data_dir, 'random_clifford')
vqe_data_dir = os.path.join(cwd, 'vqe_data')
basic_vqe_data_dir = os.path.join(vqe_data_dir, 'basic')
unitary_partitioning_vqe_data_dir = os.path.join(vqe_data_dir, 'unitary_partitioning')
exact_data_dir = os.path.join(cwd, 'exact_data')
random_pauli_10000_dir = os.path.join(shadow_data_dir, "pauli_10000")
random_clifford_10000_dir = os.path.join(shadow_data_dir, "clifford_10000")

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

total_measurement_budget = 1000

for filename in ham_list[25:30]:

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

    # # Check that gs_psi_tap is the true gs
    # nrg_calc = 0
    # for p in range(H_taper.n_terms):
    #     nrg_calc += H_taper.coeff_vec[p] * single_term_expval(H_taper[p], gs_psi_tap)
    # print(nrg_calc)

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

    # Run a quantum experiment to estimate gs energy given gs_psi_tap and H_taper
    sp_circ_copy = sp_circ.copy()
    gs_nrg_basic_vqe = measureObservable(sp_circ_copy, range(sp_circ_copy.num_qubits), H_taper, shots=total_measurement_budget, total_shots_provided=True)
    print("gs_nrg_basic_vqe = " , gs_nrg_basic_vqe)

    with open(os.path.join(basic_vqe_data_dir, filename), 'w') as file:
        data = {"total measurement budget" : total_measurement_budget, 
                "gs_nrg_tap" : str(gs_nrg_tap),
                "gs_nrg_basic_vqe" : str(gs_nrg_basic_vqe)}
        json.dump(data,file)

    # Run a quantum experiment to estimate gs energy given gs_psi_tap and H_taper using unitary partitioning
    clique_cover = H_taper.clique_cover(edge_relation='AC')
    gs_nrg_up_vqe = 0
    for clique in clique_cover.values():
        ac_op = AntiCommutingOp.from_PauliwordOp(clique)
        pop,rots,w,op = ac_op.unitary_partitioning()
        op_mat = op.to_sparse_matrix
        sp_circ_copy = sp_circ.copy()
        if rots != None:
            rots = rots[::-1]
            for rot, coeff in rots:
                rot = rot * (-0.5*coeff)
                rot = PauliwordOp_to_QuantumCircuit(rot)           
                sp_circ_copy = sp_circ_copy.compose(rot, list(range(sp_circ_copy.num_qubits)))
        gs_nrg_up_vqe += w * measureObservable(sp_circ_copy, range(sp_circ_copy.num_qubits), pop, shots=total_measurement_budget, total_shots_provided=True)
    print("gs_nrg_up_vqe = " , gs_nrg_up_vqe)
    with open(os.path.join(unitary_partitioning_vqe_data_dir, filename), 'w') as file:
        data = {"total_measurement_budget" : total_measurement_budget,
                "gs_nrg_tap" : str(gs_nrg_tap),
                "gs_nrg_up_vqe" : str(gs_nrg_up_vqe)}
        json.dump(data, file)

    # Estimate gs energy given gs_psi_tap and H_taper using classical shadows
    sp_circ_copy = sp_circ.copy()
    H_taper_dict = H_taper.to_dictionary
    obs = list(H_taper_dict.keys())
    classical_shadow = ClassicalShadow(sp_circ_copy, obs)
    classical_shadow.createClassicalShadows(unitary_ensemble="random clifford", num_shadows=total_measurement_budget)

    classical_shadow.observables = [PauliwordOp.from_dictionary({o : w}).to_sparse_matrix for o,w in zip(list(H_taper_dict.keys()), list(H_taper_dict.values()))]
    obs, results = classical_shadow.linearPredictions(10)
    gs_nrg_tap_cst = 0
    for w, exp in zip(list(H_taper_dict.values()), results):
        gs_nrg_tap_cst += exp 

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

    with open(os.path.join(random_clifford_10000_dir, filename), 'w') as file:
        shadow_data = {}
        for i in range(len(classical_shadow.shadows)):
            shadow = classical_shadow.shadows[i]
            # shadow_data[i] = {"data": str(shadow.data.tolist()), "indices": str(shadow.nonzero()[0].tolist()), "indptr": str(shadow.nonzero()[1].tolist())}
        data = {"num_shadows" : total_measurement_budget,
                "gs_nrg_tap" : str(gs_nrg_tap),
                "gs_nrg_tap_cst" : str(gs_nrg_tap_cst),
                "gs_nrg_tap_cst_up" : str(gs_nrg_tap_cst_up)
                }
        json.dump(data, file)

    # Estimate gs energy given gs_nrg_tap and H_taper using classical shadows in the Pauli basis
    sp_circ_copy = sp_circ.copy()
    H_taper_dict = H_taper.to_dictionary
    obs = list(H_taper_dict.keys())
    classical_shadow = ClassicalShadow(sp_circ_copy, obs)
    classical_shadow.createClassicalShadows(unitary_ensemble="pauli", num_shadows=1000)

    classical_shadow.observables = [PauliwordOp.from_dictionary({o : w}).to_sparse_matrix for o,w in zip(list(H_taper_dict.keys()), list(H_taper_dict.values()))]
    obs, results = classical_shadow.linearPredictions(10)
    gs_nrg_tap_cst = 0
    for w, exp in zip(list(H_taper_dict.values()), results):
        gs_nrg_tap_cst += exp 

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

    with open(os.path.join(random_pauli_10000_dir, filename), 'w') as file:
        shadow_data = {}
        for i in range(len(classical_shadow.shadows)):
            shadow = classical_shadow.shadows[i]
            # shadow_data[i] = {"data": str(shadow.data.tolist()), "indices": str(shadow.nonzero()[0].tolist()), "indptr": str(shadow.nonzero()[1].tolist())}
        data = {"num_shadows" : total_measurement_budget,
                "gs_nrg_tap" : str(gs_nrg_tap),
                "gs_nrg_tap_cst" : str(gs_nrg_tap_cst),
                "gs_nrg_tap_cst_up" : str(gs_nrg_tap_cst_up)
                }
        json.dump(data, file)

#### N.B. to reconstruct the csr_matrix for each shadow:
#### mat = csr_matrix((obj['data'], (obj['indices'], obj['indptr'])))