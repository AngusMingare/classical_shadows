from qiskit import QuantumCircuit, ClassicalRegister, execute, Aer
import numpy as np

def quantum_state_tomography(prep_circ, qubit, shots=1024):
    """Perform quantum state tomography on the state prepared by a quantum circuit

    Given a quantum circuit that prepares some quantum state |a>, this function
    performs quantum state tomography on a given qubit to construct its density matrix

    Parameters
    -----------
    prep_circ : QuantumCircuit
        A Qiskit QuantumCircuit object that prepares the state of interest |a> on any number of qubits
    qubit : int
        The index of the qubit to perform tomography on
    shots : int
        The number of shots used per circuit. Default is 1024

    Returns
    --------
    state_dm : np.ndarray
        The density matrix of the given qubit after the given ciruict
    """
    backend = Aer.get_backend("qasm_simulator")
    Id = np.array([[1,0],[0,1]])
    X = np.array([[0,1],[1,0]])
    Y = np.array([[0,-1j],[1j,0]])
    Z = np.array([[1,0],[0,-1]])

    # X 
    x_circ = QuantumCircuit(prep_circ.num_qubits)
    x_circ.append(prep_circ, list(range(prep_circ.num_qubits)))
    x_circ.h(qubit)
    
    x_circ.add_register(ClassicalRegister(1))
    x_circ.measure(qubit, 0)
    
    jobx = execute(x_circ, backend=backend, shots=shots)
    resultx = jobx.result()
    countsx = resultx.get_counts()

    countsx0 = countsx["0"] if "0" in countsx else 0
    countsx1 = countsx["1"] if "1" in countsx else 0
    exp_X = countsx0/shots - countsx1/shots
    
    # Y
    y_circ = QuantumCircuit(prep_circ.num_qubits)
    y_circ.append(prep_circ, list(range(prep_circ.num_qubits)))
    y_circ.sdg(qubit)
    y_circ.h(qubit)
    
    y_circ.add_register(ClassicalRegister(1))
    y_circ.measure(qubit, 0)
    
    joby = execute(y_circ, backend=backend, shots=shots)
    resulty = joby.result()
    countsy = resulty.get_counts()
    
    countsy0 = countsy["0"] if "0" in countsy else 0
    countsy1 = countsy["1"] if "1" in countsy else 0
    exp_Y = countsy0/shots - countsy1/shots
    
    
    # Z
    z_circ = QuantumCircuit(prep_circ.num_qubits)
    z_circ.append(prep_circ, list(range(prep_circ.num_qubits)))
    
    z_circ.add_register(ClassicalRegister(1))
    z_circ.measure(qubit, 0)
    
    jobz = execute(z_circ, backend=backend, shots=shots)
    resultz = jobz.result()
    countsz = resultz.get_counts()
    
    countsz0 = countsz["0"] if "0" in countsz else 0
    countsz1 = countsz["1"] if "1" in countsz else 0
    exp_Z = countsz0/shots - countsz1/shots
    
    # Reconstruct
    state_dm = Id/2 + exp_X*X/2 + exp_Y*Y/2 + exp_Z*Z/2
    return state_dm