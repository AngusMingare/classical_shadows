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
import matplotlib.pyplot as plt

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
for ham_idx in range(len(ham_list)):
    ham = ham_list[ham_idx]
    print(ham_idx, ham)

num_shots = [1000,3000,5000,7000,9000]
random_energy = [-74.2228698179132, -99.51009680555785, -99.42636149359726, -102.8588548172837, -95.5928116881379]
biased_energy = [-98.63015854953343, -98.4376864437716, -98.54725922358358, -98.62320486316898, -98.56038174007976]
fci = -98.6033017772588
hf = -98.57101106797646

# for total_measurement_budget in [1000,3000,5000,7000,9000]:

#     # total_measurement_budget = 10000

# # for filename_idx in range(40):

#     filename = ham_list[36]

#     if "SINGLET" in filename:
#         pass
#     else:
#         continue

#     if "STO-3G" in filename:
#         pass
#     else:
#         continue

#     with open(os.path.join(ham_data_dir, filename), 'r') as infile:
#         data_dict = json.load(infile)
#     H = PauliwordOp.from_dictionary(data_dict['hamiltonian'])

#     # Tapering
#     QT = QubitTapering(H)
#     hf_state   = QuantumState(np.asarray(data_dict['data']['hf_array'])) # Hartree-Fock state
#     hf_energy  = data_dict['data']['calculated_properties']['HF']['energy']
#     fci_energy = data_dict["data"]["calculated_properties"]["FCI"]["energy"]
#     H_taper   = QT.taper_it(ref_state=hf_state) 
#     print(H_taper)

#     # Exact Energy
#     gs_nrg_tap, gs_psi_tap = exact_gs_energy(H_taper.to_sparse_matrix)
#     #gs_nrg_tap = H_taper.expval(gs_psi_tap)

#     # Get qiskit circuit preparing gs_psi_tap
#     sp_circ = QuantumCircuit(gs_psi_tap.n_qubits)
#     gs_psi_tap_dict = gs_psi_tap.to_dictionary
#     gs_psi_tap_array = np.array([0]*(2**gs_psi_tap.n_qubits), dtype=complex)
#     def bin_to_int(b):
#         l = len(b)
#         val = 0
#         for a in range(l):
#             val += int(b[a]) * (2**(l-1-a))
#         return val
#     for p, w in zip(list(gs_psi_tap_dict.keys()), list(gs_psi_tap_dict.values())):
#         val = bin_to_int(p)
#         gs_psi_tap_array[val] = w 
#     sp_circ.prepare_state(gs_psi_tap_array)
#     gs_psi_tap_array = csr_array(gs_psi_tap_array).transpose()

#     clique_cover = H_taper.clique_cover(edge_relation='AC')
#     commuting_clique_cover = H_taper.clique_cover(edge_relation='C')


#     # # Estimate gs energy given gs_psi_tap and H_taper using classical shadows
#     sp_circ_copy = sp_circ.copy()
#     H_taper_dict = H_taper.to_dictionary
#     # print("Hamiltonian:")
#     # print(H_taper)
#     obs = list(H_taper_dict.keys())
    
#     # sp_circ_copy = QuantumCircuit(2)
#     # sp_circ_copy.h(0)
#     # sp_circ_copy.h(1)
#     # obs = ["XX", "IZ", "XI", "XZ"]
#     # H_taper = PauliwordOp.from_dictionary({"XX" : 0.5, "IZ" : 0.3, "XI" : 0.2, "XZ" : 0.1})
#     # H_taper_dict = H_taper.to_dictionary
#     classical_shadow = ClassicalShadow(sp_circ_copy, obs, H_taper)
#     classical_shadow.createClassicalShadows(unitary_ensemble="biased clifford", num_shadows=total_measurement_budget)

#     gs_nrg_tap_cst = 0

#     # shadows_bitstrings = [classical_shadow.shadows[i][1] for i in range(total_measurement_budget)]
#     stabilizer_rotations_list = classical_shadow.biased_clifford_stabilizer_rotations 
#     # stabilizers_bitstrings = [(stabilizer_rotations_list[i], shadows_bitstrings[i]) for i in range(len(shadows_bitstrings))]
#     # terms = [PauliwordOp.from_dictionary({list(H_taper_dict.keys())[i] : list(H_taper_dict.values())[i]}) for i in range(len(list(H_taper_dict.keys())))]

#     shadows_bitstrings = [classical_shadow.shadows[i][1] for i in range(total_measurement_budget)]
#     clique_list = classical_shadow.biased_clifford_cliques
#     cliques_bitstrings = [(clique_list[i], shadows_bitstrings[i]) for i in range(len(shadows_bitstrings))]
#     generators = [classical_shadow.biased_clifford_generators[i] for i in range(len(shadows_bitstrings))]
#     stabilizers_generatorss_bitstrings = [(stabilizer_rotations_list[i], generators[i], shadows_bitstrings[i]) for i in range(len(shadows_bitstrings))]
#     stabilizers_bitstrings = [(stabilizer_rotations_list[i], shadows_bitstrings[i]) for i in range(len(shadows_bitstrings))]
#     terms = [PauliwordOp.from_dictionary({list(H_taper_dict.keys())[i] : list(H_taper_dict.values())[i]}) for i in range(len(list(H_taper_dict.keys())))]

#     def parity(bitstring, sign_flip=False):
#         """
#         Returns 1 if the bitstring has an even number of 1s and returns 0 otherwise
#         """
#         parity = 1
#         for bit in bitstring:
#             if bit == "1":
#                 parity = parity * (-1)
#         if sign_flip:
#             parity = parity * (-1)
#         return parity

#     for term in terms:
#         # print("term")
#         # print(term)
#         pauli_string = list(term.to_dictionary.keys())[0]
#         weight = list(term.to_dictionary.values())[0]
#         shots = 0
#         counts0 = 0
#         counts1 = 0
#         for stabilizer_rotations, generators, bitstring in stabilizers_generatorss_bitstrings:
#             # print("rotation")
#             # print(stabilizer_rotations)
#             rotated_term_dict = term.perform_rotations(stabilizer_rotations).to_dictionary
#             rotated_term = list(rotated_term_dict.keys())[0]
#             rotated_weight = list(rotated_term_dict.values())[0]
#             # print(generators.perform_rotations(stabilizer_rotations))

#             flip_sign = False if rotated_weight == weight else True

#             recon_matrix, successful_mask = term.generator_reconstruction(generators)

#             Z_idxs = []
#             for p_idx in range(len(rotated_term)):
#                 p = rotated_term[p_idx]
#                 if p == "Z":
#                     Z_idxs.append(p_idx)

#             if successful_mask[0]:
#                 shots += 1
#                 subbitstring = ""
#                 for idx in Z_idxs:
#                     subbitstring += bitstring[idx]
#                 if parity(subbitstring, flip_sign) == 1:
#                     counts1 += 1
#                 else:
#                     counts0 += 1

#         if shots > 0:
#             gs_nrg_tap_cst += (counts1/shots - counts0/shots) * weight

#     classical_shadow = ClassicalShadow(sp_circ_copy, obs, H_taper)
#     classical_shadow.createClassicalShadows(unitary_ensemble="random clifford", num_shadows=total_measurement_budget)
#     classical_shadow.observables = [PauliwordOp.from_dictionary({o : w}).to_sparse_matrix for o,w in zip(list(H_taper_dict.keys()), list(H_taper_dict.values()))]
#     obs, results = classical_shadow.linearPredictions(1)
#     gs_nrg_tap_cst_rand = 0
#     for w, exp in zip(list(H_taper_dict.values()), results):
#         gs_nrg_tap_cst_rand += exp 

#     num_shots.append(total_measurement_budget)
#     random_energy.append(gs_nrg_tap_cst_rand)
#     biased_energy.append(gs_nrg_tap_cst)
#     hf = hf_energy
#     fci = fci_energy

#     # print(filename_idx)
#     print(filename)
#     print("total measurement budget = ", total_measurement_budget)
#     print("hf energy = ", hf_energy)
#     print("fci energy = ", fci_energy)
#     print("gs_nrg_tap = ", gs_nrg_tap)
#     print("gs_nrg_tap_cst = ", str(gs_nrg_tap_cst))
#     print("gs_nrg_tap_cst_rand = ", str(gs_nrg_tap_cst_rand))

#     # # Estimate gs energy given gs_psi_tap and H_taper using unitary partitioning and classical shadows
#     # classical_shadow.observables = []
#     # weights = []
#     # for clique in clique_cover.values():
#     #     ac_op = AntiCommutingOp.from_PauliwordOp(clique)
#     #     pop,rots,w,op = ac_op.unitary_partitioning()
#     #     classical_shadow.observables.append(op.to_sparse_matrix)
#     #     weights.append(w)
#     # _, results = classical_shadow.linearPredictions(10)
#     # gs_nrg_tap_cst_up = 0
#     # for w, exp in zip(weights, results):
#     #     gs_nrg_tap_cst_up += w * exp 
#     # print("gs_nrg_tap_cst_up = ", gs_nrg_tap_cst_up)

#     # # Estimate gs energy given gs_psi_tap and H_taper using commuting grouping and classical shadows
#     # classical_shadow.observables = []
#     # for clique in commuting_clique_cover.values():
#     #     classical_shadow.observables.append(clique.to_sparse_matrix)
#     # _, res = classical_shadow.linearPredictions(10)
#     # gs_nrg_tap_cst_gc = str(sum(res))
#     # print("gs_nrg_tap_cst_gc = ", gs_nrg_tap_cst_gc)

#     # with open(os.path.join(biased_clifford_data_dir, filename), 'w') as file:
#     #     shadow_data = {}
#     #     for i in range(len(classical_shadow.shadows)):
#     #         matrix, bitstring = classical_shadow.shadows[i]
#     #         matrix = matrix.qasm()
#     #         shadow_data[i] = (matrix, bitstring)
#     #     data = {"num_shadows" : total_measurement_budget,
#     #             "gs_nrg_tap" : str(gs_nrg_tap),
#     #             "gs_nrg_tap_cst" : str(gs_nrg_tap_cst),
#     #             "gs_nrg_tap_cst_rand" : str(gs_nrg_tap_cst_rand),
#     #             "shadow" : shadow_data
#     #             }
#     #     json.dump(data, file)

#     # classical_shadow.total_measurement_budget = 1000
#     # estimate = classical_shadow.OGMPauliEstimator("importance sampling")
#     # print("estimate = ", str(estimate))

fig, ax = plt.subplots(1,1)
ax.plot(num_shots, random_energy, c="r", marker="x")
ax.plot(num_shots, biased_energy, c="g", marker="o")
# ax.axhline(hf)
ax.axhline(fci)
ax.set_xlabel("No. Shots")
ax.set_ylabel("Energy")
ax.grid(True)
fig.set_figwidth(10.0)
fig.set_figheight(8.0)
ax.legend(['Random Clifford', 'Biased Clifford'])
plt.savefig("HF_biased_vs_random")
