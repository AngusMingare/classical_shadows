import os
import json
import numpy as np
from symmer import PauliwordOp, QuantumState, QubitTapering
from symmer.utils import exact_gs_energy
from symmer.evolution.exponentiation import *
from symmer.evolution.decomposition import *
from scipy.sparse import csr_matrix, kron
from functools import reduce

# Set up directories
cwd = os.getcwd()
file_dir =  os.path.dirname(cwd)
symmer_dir = os.path.join(cwd, 'src/symmer')
test_dir = os.path.join(symmer_dir, 'tests')
ham_data_dir = os.path.join(test_dir, 'hamiltonian_data')
vqe_data_dir = os.path.join(cwd, 'vqe_data')
basic_vqe_data_dir = os.path.join(vqe_data_dir, 'shots_10000')

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

for filename in ham_list:

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
    initial_state = gs_psi_tap.to_sparse_matrix

    h = csr_matrix([[1/np.sqrt(2), 1/np.sqrt(2)], [1/np.sqrt(2), -1/np.sqrt(2)]])
    sdg = csr_matrix([[1, 0], [0, -1j]])

    pX = h
    pY = h @ sdg
    pZ = csr_matrix([[1, 0], [0, 1]])
    
    gs_nrg_basic_vqe = 0
    for pauli, weight in H_taper.to_dictionary.items():
        to_be_tensored = [] 
        non_id = []
        for p_idx in range(len(pauli)):
            p = pauli[p_idx]
            if p == "X": 
                to_be_tensored.append(pX)
                non_id.append(p_idx)
            elif p == "Y": 
                to_be_tensored.append(pY)
                non_id.append(p_idx)
            elif p == "Z":
                to_be_tensored.append(pZ)
                non_id.append(p_idx)
            else: 
                to_be_tensored.append(pZ)
        change_of_basis = reduce(kron, to_be_tensored)
        output_state = change_of_basis @ initial_state
        num_shots = round(total_measurement_budget / len(list(H_taper.to_dictionary.keys())))
        if num_shots == 0:
            num_shots = 1
        sample_probs = [np.abs(i)**2 for i in output_state.toarray().flatten()]
        shots = np.random.choice(range(2**len(pauli)), size = num_shots, p=sample_probs)
        shots_bin = [bin(shot)[2:].zfill(len(pauli)) for shot in shots]
        reduced_shots_bin = []
        for shot_bin in shots_bin:
            reduced_shot_bin = ""
            for idx in non_id:
                reduced_shot_bin += shot_bin[idx]
            reduced_shots_bin.append(reduced_shot_bin)

        def parity(bitstring):
            """
            Returns 1 if the bitstring has an even number of 1s and returns -1 otherwise
            """
            parity = 1
            for bit in bitstring:
                if bit == "1":
                    parity = parity * (-1)
            return parity

        counts0 = 0
        counts1 = 0
        for reduced_shot_bin in reduced_shots_bin:
            if parity(reduced_shot_bin) == 1:
                counts1 += 1
            else:
                counts0 += 1

        exp = counts1/num_shots - counts0/num_shots
        gs_nrg_basic_vqe += weight * exp
    print(gs_nrg_basic_vqe)
    
    with open(os.path.join(basic_vqe_data_dir, filename), 'w') as file:
        data = {"total measurement budget" : total_measurement_budget, 
                "gs_nrg_tap" : str(gs_nrg_tap),
                "gs_nrg_basic_vqe" : str(gs_nrg_basic_vqe)}
        json.dump(data,file)