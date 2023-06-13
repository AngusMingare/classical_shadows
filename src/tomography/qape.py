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
    circ.h(0)
    myGate = UnitaryGate(unitary)
    c_myGate = myGate.control(1)
    circ.append(c_myGate, list(range(circ.num_qubits)))

    chi = quantumProcessTomography(circ, 0)
    choi = chiToChoi(chi)
    ev = choiToEvolution(choi)

    qft_tensors = qft_mpo_n_qubits(num_precision_bits, 8)
    
    initial_state = np.array([1,0,0,0])
    new_tensors = []
    for i in range(num_precision_bits):
        temp_state = np.asmatrix(initial_state.reshape((4,1)))
        for j in range(2**i):
            temp_state = np.matmul(ev, temp_state)
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
one = np.array([0,1]).reshape((2))

val = np.einsum('abc,adef,dgh,b,c,e,f,g,h', tensor1, tensor2, tensor3, zero, zero, zero, zero, one, one)
print(np.abs(val)**2)