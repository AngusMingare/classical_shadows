from qpt import *
from ..utils.qft_mpo import *
import numpy as np
from qiskit import QuantumCircuit
from qiskit.extensions import *

def quantumAssistedPhaseEstimation(unitary, num_precision_bits):
    myGate = UnitaryGate(unitary)
    c_myGate = myGate.control(1)

    circ = QuantumCircuit(1+int(np.log2(unitary.shape[0])))
    circ.append(c_myGate)

    chi = quantumProcessTomography(circ, 0)
    choi = chiToChoi(chi)
    ev = choiToEvolution(choi)

    qft_tensors = qft_mpo.qft_mpo_n_qubits(num_precision_bits, 8)
    
    initial_state = np.array([1,0,0,0])
    for i in range(num_precision_bits):
        temp_state = np.asmatrix(initial_state.reshape((4,1)))
        for j in range(2**i):
            temp_state = np.matmul(ev, temp_state)
        temp_tensor = qft_tensors[i]
        if i == 0:
            qft_tensors[i] = np.einsum('ab,cda->cdb', temp_state, temp_tensor)
        elif i == num_precision_bits-1:
            qft_tensors[i] = np.einsum('ab,cda->cdb', temp_state, temp_tensor)
        else:
            qft_tensors[i] = np.einsum('ab,cdea->cdeb', temp_state, temp_tensor)
        
    return qft_tensors

u = np.array([[1,0],[0,np.exp(2*np.pi*1j*(1/8))]])
tensors = quantumAssistedPhaseEstimation(u, 3)

tensor1 = tensors[0]
tensor2 = tensors[1]
tensor3 = tensors[2]

zero = np.array([1,0]).reshape((2))
one = np.array([0,1]).reshape((2))

val = np.einsum('abc,adef,dgh,b,c,e,f,g,h', tensor1, tensor2, tensor3, zero, zero, zero, zero, one, one)
print(val)