import sys
sys.path.append('../utils')
from qft_mpo import *
from qpt import *
import numpy as np
from qiskit import QuantumCircuit
from qiskit.extensions import *

def quantumAssistedPhaseEstimation(unitary, state_prep_circ, num_precision_bits):
    circ = QuantumCircuit(1+int(np.log2(unitary.shape[0])))
    circ.append(state_prep_circ, list(range(1, circ.num_qubits)))
    myGate = UnitaryGate(unitary)
    c_myGate = myGate.control(1)
    circ.append(c_myGate, list(range(circ.num_qubits)))

    chi = quantumProcessTomography(circ, 0)
    choi = chiToChoi(chi)
    ev = choiToEvolution(choi)
    print(ev)

    qft_tensors = qft_mpo_n_qubits(num_precision_bits, 8)
    
    initial_state = np.array([1/2,1/2,1/2,1/2])
    new_tensors = []
    for i in range(num_precision_bits):
        temp_state = np.asmatrix(initial_state.reshape((4,1)))
        for _ in range(2**i):
            temp_state = np.matmul(ev, temp_state)
        print(temp_state) 
        temp_state = np.asmatrix(np.array([[temp_state[0,0], temp_state[2,0]],[temp_state[1,0], temp_state[3,0]]]))
        temp_tensor = qft_tensors[i]
        if i == 0:
            new_tensors.append(np.einsum('ab,cda->cdb', temp_state, temp_tensor))
        elif i == num_precision_bits-1:
            new_tensors.append(np.einsum('ab,cda->cdb', temp_state, temp_tensor))
        else:
            new_tensors.append(np.einsum('ab,cdea->cdeb', temp_state, temp_tensor))
        
    return new_tensors

u = np.array([[1,0],[0,np.exp(2*np.pi*1j*(1/8))]])
state_prep_circ = QuantumCircuit(1)
state_prep_circ.x(0)
tensors = quantumAssistedPhaseEstimation(u, state_prep_circ, 3)

tensor1 = tensors[0]
tensor2 = tensors[1]
tensor3 = tensors[2]

zero = np.array([1,0]).reshape((2))

output_statevec = np.einsum('abc,adef,dgh,c,f,h->beg', tensor1, tensor2, tensor3, zero, zero, zero)
outcome = np.reshape(output_statevec, (8))
print("000 = ", np.abs(outcome[0])**2)
print("001 = ", np.abs(outcome[1])**2)
print("010 = ", np.abs(outcome[2])**2)
print("011 = ", np.abs(outcome[3])**2)
print("100 = ", np.abs(outcome[4])**2)
print("101 = ", np.abs(outcome[5])**2)
print("110 = ", np.abs(outcome[6])**2)
print("111 = ", np.abs(outcome[7])**2)