import os
from symmer import PauliwordOp, QuantumState, QubitTapering
from symmer.operators import AntiCommutingOp, single_term_expval
from symmer.utils import exact_gs_energy
from symmer.evolution.exponentiation import *
import json
import numpy as np
from vqe.measurements import measureObservable
from tomography.cst import ClassicalShadow
from qiskit import QuantumCircuit
from qiskit.extensions import UnitaryGate

# Set up directories
cwd = os.getcwd()
file_dir =  os.path.dirname(cwd)
symmer_dir = os.path.join(cwd, 'symmer')
test_dir = os.path.join(symmer_dir, 'tests')
ham_data_dir = os.path.join(test_dir, 'hamiltonian_data')

# Get a Hamiltonian
ham_list = os.listdir(ham_data_dir)
filename = ham_list[7]
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
print("The exact gs energy is " + str(gs_nrg_tap))

# Check that gs_psi_tap is the true gs
# nrg_calc = 0
# for p in range(H_taper.n_terms):
#     nrg_calc += H_taper.coeff_vec[p] * single_term_expval(H_taper[p], gs_psi_tap)
# print(nrg_calc)

# Get qiskit circuit preparing gs_psi_tap
sp_circ = QuantumCircuit(gs_psi_tap.n_qubits)
gs_psi_tap_dict = gs_psi_tap.to_dictionary
gs_psi_tap_array = np.array([0]*(2**gs_psi_tap.n_qubits), dtype=complex)
def bin_to_int(b):
    b = b[::-1]
    l = len(b)
    val = 0
    for a in range(l):
        val += int(b[a]) * (2**(l-1-a))
    return val
for p, w in zip(list(gs_psi_tap_dict.keys()), list(gs_psi_tap_dict.values())):
    val = bin_to_int(p)
    gs_psi_tap_array[val] = w 
sp_circ.prepare_state(gs_psi_tap_array)

# Run a quantum experiment to estimate gs energy given gs_psi_tap and H_taper
# sp_circ_copy = sp_circ.copy()
# gs_nrg_basic_vqe = measureObservable(sp_circ_copy, range(sp_circ_copy.num_qubits), H_taper, shots=1000)
# print("basic VQE: " + str(gs_nrg_basic_vqe))

# Run a quantum experiment to estimate gs energy given gs_psi_tap and H_taper using unitary partitioning
# clique_cover = H_taper.clique_cover(edge_relation='AC')
# gs_nrg_up_vqe = 0
# for clique in clique_cover.values():
#     ac_op = AntiCommutingOp.from_PauliwordOp(clique)
#     pop,rots,w,op = ac_op.unitary_partitioning()
#     sp_circ_copy = sp_circ.copy()
#     if rots != None:
#           for rot, coeff in rots:
#             rot = rot * coeff
#             rot = rot.to_qiskit
#             exp_rot_mat = rot.exp_i()
#             sp_circ_copy.append(exp_rot_mat, list(range(sp_circ_copy.num_qubits)))
#     gs_nrg_up_vqe += w * measureObservable(sp_circ_copy, range(sp_circ_copy.num_qubits), pop, shots=1000)
# print("UP VQE: " + str(gs_nrg_up_vqe))

# Estimate gs energy given gs_psi_tap and H_taper using classical shadows
sp_circ_copy = sp_circ.copy()
H_taper_dict = H_taper.to_dictionary 
obs = list(H_taper_dict.keys())
classical_shadow = ClassicalShadow(sp_circ_copy, obs)
classical_shadow.createClassicalShadows(unitary_ensemble="random clifford", num_shadows=1000)
obs, results = classical_shadow.linearPredictions(10)

for o, exp in (obs, results):
    print(o)
    print(exp)

# Estimate gs energy given gs_psi_tap and H_taper using unitary partitioning and classical shadows
