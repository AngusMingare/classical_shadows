from qiskit import Aer, ClassicalRegister, execute
from symmer import PauliwordOp
import numpy as np
from itertools import product

def measurePauliString(circ, qubit_indices, pauli_string):
    """Add circuit measurements given by the Pauli string
    
    This function appends some single-qubit gates to the circ so that
    the circuit measures in the bases given by the pauli_string on 
    the given qubit_indices.

    Parameters
    -----------
    circ : QuantumCircuit
        The QISKIT quantum circuit object
    qubit_indices : list of int
        The indices of the qubits associated with the Pauli string measurement
    pauli_string : str
        String composed of I,X,Y,Z corresponding to the observable to be measured
    
    Returns
    --------
    circ_out : circuit
        The original circ updated with the Pauli string basis measurement
    """

    assert(len(qubit_indices) == len(pauli_string))

    qubits_to_measure = []
    cl_bits = list(range(circ.num_clbits))

    for idx, p in zip(qubit_indices, list(pauli_string.strip())):
        idx = len(pauli_string)-1-idx
        if p == "I":
            pass
        elif p == "X":
            circ.h(idx)
            qubits_to_measure.append(idx)
        elif p == "Y":
            circ.sdg(idx)
            circ.h(idx)
            qubits_to_measure.append(idx)
        elif p == "Z":
            qubits_to_measure.append(idx)
    
    diff = len(qubits_to_measure) - len(cl_bits)
    if diff <= 0:
        circ.measure(qubits_to_measure, cl_bits[:len(qubits_to_measure)])
    else:
        cr = ClassicalRegister(diff)
        circ.add_register(cr)
        circ.measure(qubits_to_measure, range(circ.num_clbits))
    
    return circ

def measurePauliObservable(circ, qubit_indices, pauli_observable, shots=1024):
    """Measure an observable given by a single Pauli string
    
    This function calculates the expectation value <psi | O | psi>
    where |psi> is encoded by the input circ and O is given by a 
    single Pauli string (pauli_observable)

    Parameters
    -----------
    circ : QuantumCircuit
        The QISKIT quantum circuit object
    qubit_indices : list of int
        The indices of the qubits associated with the Pauli string measurement
    pauli_observable : str
        String composed of I,X,Y,Z corresponding to the observable to be measured
    shots : int
        Number of shots for each quantum circuit
    
    Returns
    --------
    expectation_value : float
        The value <psi | pauli_observable | psi>
    """

    circ = measurePauliString(circ, qubit_indices=qubit_indices, pauli_string=pauli_observable)

    backend = Aer.get_backend("qasm_simulator")
    job = execute(circ, backend=backend, shots=shots)
    result = job.result()
    counts = result.get_counts()

    def parity(bitstring):
        """
        Returns 1 if the bitstring has an even number of 1s and returns 0 otherwise
        """
        parity = 1
        for bit in bitstring:
            if bit == "1":
                parity = parity * (-1)
        return parity

    counts0 = 0
    counts1 = 0
    for bitstring, freq in counts.items():
        if parity(bitstring) == 1:
            counts1 += freq
        else:
            counts0 += freq
    
    return counts1/shots - counts0/shots 

def decomposeObservable(observable):
    """Decompose an observable into the Pauli basis

    Given a Hermitian operator O we can decompose it into a sum of Pauli strings.

    Parameters
    -----------
    observable : np ndarray
        The Hermitian observable to be decomposed

    Returns
    --------
    pauliSum : dict
        A dictionary of the form {pauliString : coefficient} so that observable = sum(coefficient * pauliString)
    """
    assert(isinstance(observable, np.ndarray))
    observable = np.asmatrix(observable)
    # assert(observable.H == observable)

    # Define Pauli matrices
    sigma_x = np.array([[0, 1], [1, 0]])
    sigma_y = np.array([[0, -1j], [1j, 0]])
    sigma_z = np.array([[1, 0], [0, -1]])
    identity = np.eye(2)

    # Get the number of qubits
    n_qubits = int(np.log2(observable.shape[0]))

    # Create an empty dictionary for the Pauli basis coefficients
    pauli_dict = {}

    # Iterate over all Pauli basis matrices
    for pauli_indices in product([0, 1, 2, 3], repeat=n_qubits):
        pauli_string = ""
        pauli_matrix = np.eye(1, dtype=np.complex128)
        for index in pauli_indices:
            if index == 0:
                pauli_string += "I"
                pauli_matrix = np.kron(pauli_matrix, identity)
            elif index == 1:
                pauli_string += "X"
                pauli_matrix = np.kron(pauli_matrix, sigma_x)
            elif index == 2:
                pauli_string += "Y"
                pauli_matrix = np.kron(pauli_matrix, sigma_y)
            elif index == 3:
                pauli_string += "Z"
                pauli_matrix = np.kron(pauli_matrix, sigma_z)
        coeff = np.trace(np.matmul(observable, pauli_matrix)) / (2 ** n_qubits)
        pauli_dict[pauli_string] = coeff
    
    return pauli_dict 

def measureObservable(circ, qubit_indices, observable, shots=1024, total_shots_provided=False):
    """Measure the expectation value of an observable
    
    This function calculates the expectation value <psi | observable | psi>
    where |psi> is encoded by the input circ

    Parameters
    -----------
    circ : QuantumCircuit
        The QISKIT quantum circuit object
    qubit_indices : list of int
        The indices of the qubits associated with the Pauli string measurement
    observable : Union(np.ndarray, dict, symmer.PauliwordOp)
        Encodes the information of the observable to be measured
    shots : int
        Number of shots for each quantum circuit
    
    Returns
    --------
    expectation_value : float
        The value <psi | observable | psi>
    """

    if isinstance(observable, dict):
        pass
    elif isinstance(observable, np.ndarray):
        observable = decomposeObservable(observable)
    elif isinstance(observable, PauliwordOp):
        observable = observable.to_dictionary
    
    # <observable> = sum(coeff * <pauliString>)
    total_exp = 0
    total_num_pauli = len(list(observable.keys()))
    for pauliString, coefficient in observable.items():
        all_I = True 
        for p in pauliString:
            if p != 'I':
                all_I = False
        if all_I:
            total_num_pauli -= 1
            exp = 1
        else:
            usable_shots = round(shots / total_num_pauli) if total_shots_provided else shots
            exp = measurePauliObservable(circ=circ, qubit_indices=qubit_indices, pauli_observable=pauliString, shots=usable_shots)
        total_exp += exp*coefficient
    
    return total_exp