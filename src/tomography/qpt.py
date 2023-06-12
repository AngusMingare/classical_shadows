from qiskit import QuantumCircuit
from qst import quantumStateTomography
import numpy as np

def quantumProcessTomography(circ, qubit, shots=1024):
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
    
    qst1 = quantumStateTomography(prep1, qubit, shots)
    qst2 = quantumStateTomography(prep2, qubit, shots)
    qst3 = quantumStateTomography(prep3, qubit, shots)
    qst4 = quantumStateTomography(prep4, qubit, shots)
    
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

def chiToChoi(chi_matrix):
    """Convert chi matrix to Choi matrix

    Convert the chi matrix describing the quantum process on a single qubit
    to the corresponding Choi matrix.

    Parameters
    -----------
    chi_matrix : np.matrix
        The chi matrix of a quantum process
    
    Returns
    --------
    choi_matrix : np.matrix
        The corresponding Choi matrix for the quantum process given by chi
    """
    choi_matrix = np.zeros((4,4))
    def vectorise_2x2(mat):
        return np.array([mat[0][0], mat[1][0], mat[0][1], mat[1][1]]).reshape((4))
    
    X = np.array([[0,1],[1,0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1,0], [0,-1]])
    Id = np.array([[1,0],[0,1]])   
    paulis = [Id, X, Y, Z]
    for m in range(4):
        for n in range(4):
            pauli_m = vectorise_2x2(paulis[m])
            pauli_n = vectorise_2x2(paulis[n]).conjugate()
            to_add = chi_matrix[m,n] * np.outer(pauli_m, pauli_n)
            choi_matrix = choi_matrix + to_add
    choi_matrix = np.asmatrix(choi_matrix)

    return choi_matrix

def choiToEvolution(choi_matrix):
    """Convert a Choi matrix into an evolution matrix

    Convert a Choi matrix describing a quantum process on a single qubit
    into the corresponding evolution matrix

    Parameters
    -----------
    choi_matrix : np.matrix
        The Choi matrix of a quantum process
    
    Returns
    --------
    ev_matrix : np.matrix
        The corresponding evolution matrix for the quantum process given by Choi
    """

    ev_matrix = np.asmatrix(np.array([
        [choi_matrix[0, 0], choi_matrix[2, 0], choi_matrix[0, 2], choi_matrix[2, 2]],
        [choi_matrix[1, 0], choi_matrix[3, 0], choi_matrix[1, 2], choi_matrix[3, 2]],
        [choi_matrix[0, 1], choi_matrix[2, 1], choi_matrix[0, 3], choi_matrix[2, 3]],
        [choi_matrix[1, 1], choi_matrix[3, 1], choi_matrix[1, 3], choi_matrix[3, 3]]
    ]))

    return ev_matrix
