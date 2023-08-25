import numpy as np
from scipy.sparse import csr_matrix, kron, csr_array, linalg
from qiskit import QuantumCircuit, Aer, execute
from qiskit.extensions import UnitaryGate
from qiskit.quantum_info import random_clifford, Clifford
from functools import reduce
from typing import List
from symmer.operators import AntiCommutingOp, IndependentOp, PauliwordOp
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.opflow import I, Z, X, Y

class ClassicalShadow:
    def __init__(self, circ, observables, hamiltonian=None):
        self.allowed_measurement_bases = ["pauli", "random clifford", "biased clifford", "globally biased pauli", "locally biased pauli", "derandomised pauli"]
        self.shadows = []
        self.sparse_shadows = []
        self.num_shadows = 0 
        self.unitary_ensemble = None
        self.preparation_circuit = circ
        self.num_qubits = circ.num_qubits
        self.observables = observables
        self.hamiltonian = hamiltonian
        self.weights = []
        self.total_measurement_budget = 0

    def samplePauliEnsemble(self):
        pauli_unitary = QuantumCircuit(self.num_qubits)
        pauli_list = []
        paulis = ["x", "y", "z"]
        for i in range(self.num_qubits):
            pauli = np.random.choice(paulis)
            pauli_list.append(pauli)
            if pauli == "x":
                pauli_unitary.h(i)
            elif pauli == "y":
                pauli_unitary.sdg(i)
                pauli_unitary.h(i)
            elif pauli == "z":
                pass

        return pauli_unitary, pauli_list

    def sampleRandomCliffordEnsemble(self):
        cliff = random_clifford(self.num_qubits)
        clifford_matrix = csr_array(cliff.to_matrix())
        clifford_unitary = cliff.to_circuit()
        return clifford_unitary, clifford_matrix
    
    def sampleBiasedCliffordEnsemble(self):
        num_terms = len(list(self.hamiltonian.to_dictionary.keys()))
        num_terms_to_keep = np.random.choice(range(min(10, int(round(num_terms/4))), num_terms))
        sub_hamiltonian = {}
        kept_terms = np.random.choice(list(self.hamiltonian.to_dictionary.keys()), size=num_terms_to_keep)
        for term in kept_terms:
            sub_hamiltonian[term] = self.hamiltonian.to_dictionary[term]
        sub_hamiltonian_op = PauliwordOp.from_dictionary(sub_hamiltonian)
        clique = sub_hamiltonian_op.largest_clique('C')
        independent_op = IndependentOp.symmetry_generators(clique)
        independent_op.generate_stabilizer_rotations()
        stabilizer_rotations = independent_op.stabilizer_rotations

        # commuting_clique_cover = list(self.hamiltonian.clique_cover(edge_relation='C').values())
        # clique = np.random.choice(commuting_clique_cover)
        # independent_op = IndependentOp.symmetry_generators(clique)
        # independent_op.generate_stabilizer_rotations()
        # stabilizer_rotations = independent_op.stabilizer_rotations

        clifford_matrix = csr_array(np.eye(2**self.num_qubits))
        x_mat = csr_array(np.array([[0,1],[1,0]]))
        y_mat = csr_array(np.array([[0,-1j],[1j,0]]))
        z_mat = csr_array(np.array([[1,0],[0,-1]]))
        id_mat = csr_array(np.eye(2))
        for rot, _ in stabilizer_rotations:
            p_rot = [*list(rot.to_dictionary.keys())[0]]
            to_be_tensored = []
            for p in p_rot:
                if p == "X": to_be_tensored.append(x_mat)
                elif p == "Y": to_be_tensored.append(y_mat)
                elif p == "Z": to_be_tensored.append(z_mat)
                else: to_be_tensored.append(id_mat)
            p_rot_mat = reduce(kron, to_be_tensored)
            clifford_matrix = linalg.expm(p_rot_mat.multiply(1j * (np.pi/4))) @ clifford_matrix
        clifford_unitary = QuantumCircuit(self.num_qubits)
        unitary = UnitaryGate(clifford_matrix.todense())
        clifford_unitary.append(unitary, list(range(self.num_qubits))[::-1])

        # clifford_unitary = QuantumCircuit(self.num_qubits)
        # for rot, _ in stabilizer_rotations:
        #     p_rot = rot.to_qiskit
        #     evo = PauliEvolutionGate(p_rot, time=np.pi/4)
        #     clifford_unitary.append(evo, list(range(self.num_qubits))[::-1])
        # clifford_matrix = self.getUnitaryFromCirc(clifford_unitary)

        # clifford_unitary = QuantumCircuit(self.num_qubits)
        # for rot, _ in stabilizer_rotations:
        #     p_rot = list(rot.to_dictionary.keys())[0]
        #     qc = pauli_string_to_circ(p_rot)
        #     clifford_unitary = clifford_unitary + qc
        # clifford_matrix = self.getUnitaryFromCirc(clifford_unitary)
        
        return clifford_unitary, clifford_matrix

    def sampleGloballyBiasedPauliEnsemble(self):
        clique_cover = list(self.hamiltonian.clique_cover(edge_relation='AC').values())
        clique = np.random.choice(clique_cover)
        ac_op = AntiCommutingOp.from_PauliwordOp(clique)
        pop,_,_,_ = ac_op.unitary_partitioning()
        pauli = pop.to_qiskit
        pauli_unitary = QuantumCircuit(self.num_qubits)
        pauli_unitary.append(pauli.to_instruction(), range(self.num_qubits))
        pauli_list = [*list(pop.to_dictionary.keys())[0]]
        return pauli_unitary, pauli_list
    
    def getUnitaryFromCirc(self, circ):
        if isinstance(circ, QuantumCircuit):
            backend = Aer.get_backend("unitary_simulator")
            job = execute(circ, backend=backend)
            result = job.result()
            return csr_array(result.get_unitary(circ, 9))
        elif isinstance(circ, str):
            qc = QuantumCircuit.from_qasm_str(circ)
            backend = Aer.get_backend("unitary_simulator")
            job = execute(qc, backend=backend)
            result = job.result()
            return csr_array(result.get_unitary(qc, 18))
            
    def getSnapshots(self, unitary_ensemble, num_shadows):
        if not isinstance(unitary_ensemble, str):
            raise TypeError("unitary_ensemble: please input a valid string.")
        
        if unitary_ensemble.lower() not in self.allowed_measurement_bases:
            raise ValueError(unitary_ensemble + " is not an allowed measurement basis.")
        
        self.unitary_ensemble = unitary_ensemble
        self.num_shadows = num_shadows

        sp_unitary = self.getUnitaryFromCirc(self.preparation_circuit)
        input_state = np.zeros((2**self.num_qubits, 1))
        input_state[0] = 1
        input_state = csr_array(input_state, shape=(2**self.num_qubits, 1))
        sp_state = sp_unitary @ input_state

        snapshots = []
        for _ in range(num_shadows):
            # Select unitary from ensemble
            if unitary_ensemble == "pauli":
                unitary, pauli_list = self.samplePauliEnsemble()
                matrix = self.getUnitaryFromCirc(unitary)
            elif unitary_ensemble == "globally biased pauli":
                unitary, pauli_list = self.sampleGloballyBiasedPauliEnsemble()
                matrix = self.getUnitaryFromCirc(unitary)
            elif unitary_ensemble == "random clifford":
                unitary, matrix = self.sampleRandomCliffordEnsemble()
            elif unitary_ensemble == "biased clifford":
                idx = np.random.choice([0,1], p=[1.0, 0.0])
                if idx == 0:
                    unitary, matrix = self.sampleBiasedCliffordEnsemble()
                else:
                    unitary, matrix = self.sampleRandomCliffordEnsemble()

            def get_shot(u):
                final_state = u @ sp_state
                sample_probs = [np.abs(i)**2 for i in final_state.toarray().flatten()]
                shot = np.random.choice(range(2**self.num_qubits), p=sample_probs)
                shot_bin = bin(shot)[2:].zfill(self.num_qubits)
                return shot_bin

            bitstring = get_shot(matrix)

            if unitary_ensemble == "pauli":
                snapshots.append((pauli_list, bitstring))
            elif unitary_ensemble == "random clifford":
                snapshots.append((unitary, bitstring))
            elif unitary_ensemble == "biased clifford":
                snapshots.append((unitary, bitstring))
            elif unitary_ensemble == "globally biased pauli":
                snapshots.append((pauli_list, bitstring))

        return snapshots

    def getQuantumChannel(self, unitary_ensemble):
        if not (isinstance(unitary_ensemble, str) or isinstance(unitary_ensemble, List[QuantumCircuit])):
            raise TypeError("unitary_ensemble: please input a string, a list of QuantumCircuit objects.")
        
        if isinstance(unitary_ensemble, str) and unitary_ensemble.lower() not in self.allowed_measurement_bases:
            raise ValueError(unitary_ensemble + " is not an allowed measurement basis.")
        
        if isinstance(unitary_ensemble, str) and unitary_ensemble.lower() == "random clifford":
            def channel(x):
                return (2**self.num_qubits + 1) * x - csr_matrix(np.identity(2**self.num_qubits))
            return channel
        elif isinstance(unitary_ensemble, str) and unitary_ensemble.lower() == "biased clifford":
            def channel(x):
                return (2**self.num_qubits + 1) * x - csr_matrix(np.identity(2**self.num_qubits))
            return channel
        elif isinstance(unitary_ensemble, str) and unitary_ensemble.lower() == "pauli":
            def channel(x):
                to_be_tensored = []
                for i in range(self.num_qubits):
                    m1 = 3 * x[i] - csr_matrix(np.identity(2))
                    to_be_tensored.append(m1)
                mP = reduce(kron, to_be_tensored)
                return mP
            return channel
        elif isinstance(unitary_ensemble, str) and unitary_ensemble.lower() == "globally biased pauli":
            def channel(x):
                to_be_tensored = []
                for i in range(self.num_qubits):
                    m1 = 3 * x[i] - csr_matrix(np.identity(2))
                    to_be_tensored.append(m1)
                mP = reduce(kron, to_be_tensored)
                return mP
            return channel
        else:
            #TODO
            raise NotImplementedError("Not yet.")

    def createClassicalShadows(self, unitary_ensemble, num_shadows):
        snapshots = self.getSnapshots(unitary_ensemble, num_shadows)
        self.shadows = [snapshots[i] for i in range(num_shadows)]
        return 

    def getSparseShadows(self):
            
            zero = np.array([1,0])
            one = np.array([0,1])
            h = csr_matrix(np.array([[1/np.sqrt(2), 1/np.sqrt(2)],[1/np.sqrt(2), -1/np.sqrt(2)]]))
            sdg = csr_matrix(np.array([[1,0],[0,-1j]]))

            channel = self.getQuantumChannel(self.unitary_ensemble)

            for i in range(self.num_shadows):
                if self.unitary_ensemble == "random clifford" or self.unitary_ensemble == "biased clifford":
                    circ, bitstring = self.shadows[i]
                    matrix = self.getUnitaryFromCirc(circ)
                    def bitstring_to_array(bitstring):
                        temp = [zero if bit == '0' else one for bit in bitstring]
                        vec = reduce(np.kron, temp)
                        return csr_matrix(np.outer(vec, vec))
                    bitstring_array = bitstring_to_array(bitstring)
                    sparse_snapshot = matrix.getH() @ bitstring_array @ matrix
                    sparse_shadow = channel(sparse_snapshot)
                    self.sparse_shadows.append(sparse_shadow)

                elif self.unitary_ensemble == "pauli" or self.unitary_ensemble == "globally biased pauli":
                    pauli_list, bitstring = self.shadows[i]

                    bitstring = bitstring[::-1]
                    def get_sq_snapshots(bitstring, p_list):
                        bitvals = [csr_matrix(np.outer(zero,zero)) if bit == '0' else csr_matrix(np.outer(one,one)) for bit in bitstring]
                        def pauli_to_mat(p):
                            if p.lower() == 'x': return h
                            elif p.lower() == 'y': return h @ sdg 
                            else: return csr_matrix(np.eye(2))
                        us = [pauli_to_mat(p) for p in p_list]
                        udgs = [u.getH() for u in us]
                        snapshot = [csr_matrix(udgs[j] @ bitvals[j] @ us[j]) for j in range(self.num_qubits)]
                        return snapshot
                    sparse_snapshot = get_sq_snapshots(bitstring, pauli_list)
                    sparse_shadow = channel(sparse_snapshot)
                    self.sparse_shadows.append(sparse_shadow)
            return

    def binShadows(self, num_bins):
        self.getSparseShadows()
        size_of_bin = int(np.floor(self.num_shadows / num_bins))
        binned_shadows = []
        for i in range(num_bins-1):
            binned_shadow = (1/size_of_bin) * np.sum(self.sparse_shadows[i*size_of_bin : (i+1)*size_of_bin])
            binned_shadows.append(binned_shadow)
        remaining_shadows = self.sparse_shadows[size_of_bin * (num_bins - 1) :]
        num_remaining_shadows = len(remaining_shadows)
        binned_shadow = (1/num_remaining_shadows) * np.sum(remaining_shadows)
        binned_shadows.append(binned_shadow)
        return binned_shadows
    
    def medianOfMeansEstimate(self, shadows):
        outputs = []
        for i in range(len(self.observables)): 
            observable = self.observables[i]
            estimates = []
            for j in range(len(shadows)):
                shadow = shadows[j]
                prod = observable @ shadow
                estimate = prod.trace()
                estimates.append(estimate)
            outputs.append(np.median(estimates))

        return self.observables, outputs
    
    def linearPredictions(self, num_bins):
    
        # Construct num_bins estimators
        binned_shadows = self.binShadows(num_bins)

        # Median of means estimation
        outputs = self.medianOfMeansEstimate(binned_shadows)

        return outputs
    
    def getDerandomisedPauliMeasurements(self, measurements):
        return
    
    def getUnbiasingFunc(self, method="LBCS"):
        allowed_methods = ["LBCS", "DCS", "importance sampling"]
        if method not in allowed_methods:
            raise NotImplementedError
        
        if method == "LBCS":
            func = 1
        elif method == "DCS":
            func = 1
        elif method == "importance sampling":
            def func(chosen_pauli, pauli_string, pauli_distribution):
                chosen_pauli_prob = pauli_distribution[chosen_pauli]
                def delta(p, q):
                    if p == q:
                        return 1
                    else:
                        return 0
                
                return (1/chosen_pauli_prob) * delta(pauli_string, chosen_pauli)

        return func

    def getPauliSet(self, method="LBCS"):
        allowed_methods = ["LBCS", "DCS", "importance sampling"]
        if method not in allowed_methods:
            raise NotImplementedError 
        
        if method == "LBCS":
            return
        elif method == "DCS":
            return
        elif method == "importance sampling":
            return list(self.hamiltonian.to_dictionary.keys())
        
    def getPauliDistribution(self, method="LBCS"):
        allowed_methods = ["LBCS", "DCS", "importance sampling"]
        if method not in allowed_methods:
            raise NotImplementedError 
        
        if method == "LCBS":
            return
        elif method == "DCS":
            return
        elif method == "importance sampling":
            dist = {}
            strings = list(self.hamiltonian.to_dictionary.keys())
            weights = list(self.hamiltonian.to_dictionary.values())
            l1_norm_weights = sum([np.abs(weight) for weight in weights])
            for idx in range(len(strings)):
                dist[strings[idx]] = np.abs(weights[idx]) / l1_norm_weights
            return dist

    def OGMPauliEstimator(self, method="LBCS"):
        allowed_methods = ["LBCS", "DCS", "importance sampling"]
        if method not in allowed_methods:
            raise NotImplementedError 
        
        weights = list(self.hamiltonian.to_dictionary.values())
        pauli_strings = list(self.hamiltonian.to_dictionary.keys())

        pauli_set = self.getPauliSet(method)
        pauli_distribution = self.getPauliDistribution(method)
        f = self.getUnbiasingFunc(method)

        def supp(p):
            support = []
            for idx in range(len(p)):
                if p[idx] == "I":
                    pass
                else:
                    support.append(idx)
            return support

        zeros = np.zeros((2**self.num_qubits, 1))
        zeros[0] = 1
        zeros = csr_array(zeros)
        input_state = self.getUnitaryFromCirc(self.preparation_circuit) @ zeros

        def mu(p, supp_q):

            def get_shot(output):
                sample_probs = [np.abs(i)**2 for i in output.toarray().flatten()]
                shot = np.random.choice(range(2**self.num_qubits), p=sample_probs)
                shot_bin = bin(shot)[2:].zfill(self.num_qubits)
                return shot_bin

            def get_mu(x):
                h = csr_array(np.array([[1/np.sqrt(2), 1/np.sqrt(2)], [1/np.sqrt(2), -1/np.sqrt(2)]]))
                sdg = csr_array(np.array([[1,0],[0,-1j]]))
                measure_x = h
                measure_y = h @ sdg
                to_be_tensored = []
                for pauli_op in x:
                    if pauli_op == "X": to_be_tensored.append(measure_x)
                    elif pauli_op == "Y": to_be_tensored.append(measure_y)
                    else: to_be_tensored.append(np.eye(2))
                measurement_op = reduce(kron, to_be_tensored)
                output_state = measurement_op @ input_state
                shot = get_shot(output_state)
                return [*shot]

            m = 1
            mu_p = get_mu(p)
            for idx in range(len(supp_q)):
                mu_pidx = mu_p[idx]
                if int(mu_pidx) == 1:
                    mu_pidx = -1
                else:
                    mu_pidx = 1
                m = m * int(mu_pidx)
            
            return m
        
        measured_paulis = []
        for j in range(self.total_measurement_budget):
            chosen_pauli = np.random.choice(pauli_set, p=list(pauli_distribution.values()))
            measured_paulis.append(chosen_pauli)

        estimates = []
        for measured_pauli in measured_paulis:
            estimate = 0
            for j in range(len(weights)):
                estimate += weights[j] * f(measured_pauli, pauli_strings[j], pauli_distribution) * mu(measured_pauli, supp(pauli_strings[j]))
            estimates.append(estimate)

        final_estimate = np.mean(estimates)

        return final_estimate


def pauli_string_to_circ(pauli_string):
    pauli_string_split = [*pauli_string]
    qc = QuantumCircuit(len(pauli_string_split))
    non_id = []

    for idx in range(len(pauli_string_split)):
        pauli_op = pauli_string_split[idx]
        if pauli_op == "X":
            non_id.append(idx)
            qc.h(idx)
        elif pauli_op == "Y":
            non_id.append(idx)
            qc.sdg(idx)
            qc.h(idx)
        elif pauli_op == "Z":
            non_id.append(idx)
        else:
            pass

        if len(non_id) == 0:
            return qc
        
        for x in range(len(non_id)-1):
            idx = non_id[x]
            idxp1 = non_id[x+1]
            qc.cnot(idx, idxp1)
        
        qc.s(non_id[-1])

        for x in range(len(non_id)-1):
            len_non_id = len(non_id)
            idx = non_id[len_non_id-1-x-1]
            idxp1 = non_id[len_non_id-1-x]
            qc.cnot(idx, idxp1)

    for idx in range(len(pauli_string_split)):
        pauli_op = pauli_string_split[idx]
        if pauli_op == "X":
            qc.h(idx)
        elif pauli_op == "Y":
            qc.h(idx)
            qc.s(idx)
        else:
            pass

    return qc
        




        



if __name__ == "__main__":
    ###
    # Quick Test
    ###

    PauliX = [[0, 1], [1,0]]
    PauliZ = [[1,0], [0,-1]]

    qc = QuantumCircuit(1)
    qc.h(0)
    obs = [csr_matrix(PauliZ), csr_matrix(PauliX)]

    classical_shadow =  ClassicalShadow(qc, obs)

    classical_shadow.createClassicalShadows(unitary_ensemble="pauli", num_shadows=1000)
    obs, results = classical_shadow.linearPredictions(10)
    print(results)