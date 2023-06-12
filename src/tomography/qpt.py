from qiskit import QuantumCircuit
from qst import quantum_state_tomography
import numpy as np

def qpt(circ, qubit, shots=1024):
    """Perform quantum process tomography

    Given a quantum circuit acting on any number of qubits, determine the process
    affecting the given qubit.

    Parameters
    -----------
    circ : QuantumCircuit
        A Qiskit QuantumCircuit object
    qubit : int
        The index of the qubit to perform tomography on
    shots : int
        The number of shots used per circuit. Default is 1024

    Returns
    --------
    chi : np.matrix
        The 4x4 chi matrix describing the process on the given qubit
    """
    prep1 = QuantumCircuit(circ.num_qubits)
    prep1.append(circ, list(range(circ.num_qubits)))
    
    prep2 = QuantumCircuit(circ.num_qubits)
    prep2.x(qubit)
    prep2.append(circ, list(range(circ.num_qubits)))
    
    prep3 = QuantumCircuit(circ.num_qubits)
    prep3.h(qubit)
    prep3.append(circ, list(range(circ.num_qubits)))
    
    prep4 = QuantumCircuit(circ.num_qubits)
    prep4.h(qubit)
    prep4.s(qubit)
    prep4.append(circ, list(range(circ.num_qubits)))
    
    qst1 = quantum_state_tomography(prep1, qubit, shots)
    qst2 = quantum_state_tomography(prep2, qubit, shots)
    qst3 = quantum_state_tomography(prep3, qubit, shots)
    qst4 = quantum_state_tomography(prep4, qubit, shots)
    
    rho1 = qst1
    rho4 = qst2
    rho2 = qst3 -1j*qst4 - (1-1j)*(rho1 + rho4)/2
    rho3 = qst3 +1j*qst4 - (1+1j)*(rho1 + rho4)/2
    
    Lambda = np.asmatrix(np.array([[1/2, 0, 0, 1/2], [0, 1/2, 1/2, 0], [0, 1/2, -1/2, 0], [1/2, 0, 0, -1/2]]))
    Rhos = np.asmatrix(np.array([
        [rho1[0,0], rho1[0,1], rho2[0,0], rho2[0,1]],
        [rho1[1,0], rho1[1,1], rho2[1,0], rho2[1,1]],
        [rho3[0,0], rho3[0,1], rho4[0,0], rho4[0,1]],
        [rho3[1,0], rho3[1,1], rho4[1,0], rho4[1,1]]
    ]))
    chi = np.matmul(Lambda, np.matmul(Rhos, Lambda))
    
    return chi