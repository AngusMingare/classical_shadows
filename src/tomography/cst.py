import numpy as np
from scipy.sparse import csr_matrix, kron, csr_array
from qiskit import QuantumCircuit, Aer, execute
from qiskit.quantum_info import random_clifford
from functools import reduce
from typing import List
from symmer import PauliwordOp


class ClassicalShadow:
    def __init__(self, circ, observables):
        self.allowed_measurement_bases = ["pauli", "random clifford"]
        self.shadows = []
        self.sparse_shadows = []
        self.num_shadows = 0 
        self.unitary_ensemble = None
        self.preparation_circuit = circ
        self.num_qubits = circ.num_qubits
        self.observables = observables

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
    
    def getUnitaryFromCirc(self, circ):
        backend = Aer.get_backend("unitary_simulator")
        job = execute(circ, backend=backend)
        result = job.result()
        return csr_array(result.get_unitary(circ, 9))
    
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
            elif unitary_ensemble == "random clifford":
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
        elif isinstance(unitary_ensemble, str) and unitary_ensemble.lower() == "pauli":
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

                if self.unitary_ensemble == "random clifford":
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

                elif self.unitary_ensemble == "pauli":
                    pauli_list, bitstring = self.shadows[i]
                    circ = QuantumCircuit(self.num_qubits)
                    for i in range(len(pauli_list)):
                        pauli = pauli_list[i]
                        if pauli == "x":
                            circ.h(i)
                        elif pauli == "y":
                            circ.sdg(i)
                            circ.h(i)
                        elif pauli == "z":
                            pass
                    matrix = self.getUnitaryFromCirc(circ)
                    bitstring = bitstring[::-1]
                    def get_sq_snapshots(bitstring, p_list):
                        bitvals = [np.outer(zero,zero) if bit == '0' else np.outer(one,one) for bit in bitstring]
                        def pauli_to_mat(p):
                            if p == 'x': return h
                            elif p == 'y': return h @ sdg 
                            else: return csr_matrix(np.eye(2))
                        us = [pauli_to_mat(p) for p in p_list]
                        udgs = [u.getH() for u in us]
                        snapshot = [csr_matrix(udgs[i] @ bitvals[i] @ us[i]) for i in range(self.num_qubits)]
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