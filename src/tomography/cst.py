import numpy as np
from typing import Union, List
import symmer
from scipy.sparse import csr_matrix, kron
from qiskit import QuantumCircuit, ClassicalRegister, Aer, execute
from qiskit.quantum_info import random_clifford, Clifford
from src.vqe.measurements import decomposeObservable
from symmer import PauliwordOp


class ClassicalShadow:
    def __init__(self, circ, observables):
        self.allowed_measurement_bases = ["pauli", "random clifford"]
        self.shadows = []
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
        clifford_unitary = cliff.to_circuit()
        return clifford_unitary
    
    def getUnitaryFromCirc(self, circ):
        backend = Aer.get_backend("unitary_simulator")
        job = execute(circ, backend=backend)
        result = job.result()
        return csr_matrix(result.get_unitary(circ, 9))
    
    def getSnapshots(self, unitary_ensemble, num_shadows):
        if not isinstance(unitary_ensemble, str):
            raise TypeError("unitary_ensemble: please input a valid string.")
        
        if unitary_ensemble.lower() not in self.allowed_measurement_bases:
            raise ValueError(unitary_ensemble + " is not an allowed measurement basis.")
        
        self.unitary_ensemble = unitary_ensemble
        self.num_shadows = num_shadows

        snapshots = []
        for _ in range(num_shadows):
            # Select unitary from ensemble
            if unitary_ensemble == "pauli":
                unitary, pauli_list = self.samplePauliEnsemble()
            elif unitary_ensemble == "random clifford":
                unitary = self.sampleRandomCliffordEnsemble()

            # Make measurment
            qc = QuantumCircuit(self.num_qubits)
            qc.append(self.preparation_circuit, range(self.num_qubits))
            qc.append(unitary, range(self.num_qubits))
            cr = ClassicalRegister(self.num_qubits)
            qc.add_register(cr)
            qc.measure(range(self.num_qubits), range(self.num_qubits))
            backend = Aer.get_backend("qasm_simulator")
            job = execute(qc, backend=backend, shots=1, memory=True)
            bitstring = list(job.result().get_counts().keys())[0]

            zero = np.array([1,0])
            one = np.array([0,1])
            h = np.array([[1/np.sqrt(2), 1/np.sqrt(2)],[1/np.sqrt(2), -1/np.sqrt(2)]])
            s = np.array([[1,0],[0,1j]])
            sdg = np.array([[1,0],[0,-1j]])
 
            # Apply U^dagger
            if unitary_ensemble == "random clifford":
                unitary_inverse = unitary.inverse()
                u = self.getUnitaryFromCirc(unitary)
                udg = self.getUnitaryFromCirc(unitary_inverse) 
                temp = []
                for bit in bitstring:
                    if bit == "0":
                        temp.append(zero)
                    else:
                        temp.append(one)
                for i in range(len(temp)-1):
                    temp[i+1] = np.kron(temp[i], temp[i+1])
                vec = temp[-1]
                bitstring_array = csr_matrix(np.outer(vec, vec))
                snapshot = udg @ bitstring_array @ u
                snapshots.append(snapshot)
            elif unitary_ensemble == "pauli":
                snapshot = []
                for i in range(self.num_qubits):
                    bitstring = bitstring[::-1]
                    bit = bitstring[i]
                    sq_snapshot = np.outer(zero,zero) if bit == "0" else np.outer(one,one)
                    pauli = pauli_list[i]
                    if pauli == "x":
                        u = h
                        udg = h
                    elif pauli == "y":
                        u = h @ sdg
                        udg = s @ h
                    elif pauli == "z":
                        u = np.identity(2)
                        udg = np.identity(2)
                    sq_snapshot = udg @ sq_snapshot @ u
                    sq_snapshot = csr_matrix(sq_snapshot)
                    snapshot.append(sq_snapshot)
                snapshots.append(snapshot)
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
                for i in range(self.num_qubits-1):
                    to_be_tensored[i+1] = kron(to_be_tensored[i], to_be_tensored[i+1])
                mP = to_be_tensored[-1]
                return mP
            return channel
        else:
            #TODO
            raise NotImplementedError("Not yet.")

    def createClassicalShadows(self, unitary_ensemble, num_shadows):
        snapshots = self.getSnapshots(unitary_ensemble, num_shadows)
        print("snapshots collected")
        channel = self.getQuantumChannel(unitary_ensemble)
        self.shadows = []
        for i in range(num_shadows):
            shot = snapshots[i]
            shadow = channel(shot)
            self.shadows.append(shadow)
        print("shadows obtained")
        return 

    def binShadows(self, num_bins):
        size_of_bin = int(np.floor(self.num_shadows / num_bins))
        binned_shadows = []
        for i in range(num_bins-1):
            binned_shadow = (1/size_of_bin) * np.sum(self.shadows[i*size_of_bin : (i+1)*size_of_bin])
            binned_shadows.append(binned_shadow)
        remaining_shadows = self.shadows[size_of_bin * (num_bins - 1) :]
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