import numpy as np
from qiskit import QuantumCircuit, Aer, execute
from qiskit.visualization import plot_state_city, plot_bloch_multivector
from qiskit.ignis.verification.tomography import state_tomography_circuits, StateTomographyFitter
from qiskit.ignis.verification.tomography import process_tomography_circuits, ProcessTomographyFitter
from qiskit.quantum_info import Choi, choi2evolution

# Quantum State Tomography
def quantum_state_tomography():
    # Define the quantum circuit
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)

    # Generate the state tomography circuits
    tomo_circuits = state_tomography_circuits(qc, [0, 1])

    # Execute the tomography circuits on a simulator
    simulator = Aer.get_backend('qasm_simulator')
    job = execute(tomo_circuits, simulator, shots=10000)

    # Fit the measurement results to obtain the density matrix
    tomo_fitter = StateTomographyFitter(job.result(), tomo_circuits)
    rho = tomo_fitter.fit()

    # Print the estimated density matrix
    print(rho)

    # Plot the estimated density matrix as a cityscape
    plot_state_city(rho)


# Quantum Process Tomography
def quantum_process_tomography():
    # Define the quantum circuit
    qc = QuantumCircuit(1)
    qc.h(0)
    qc.rz(np.pi/4, 0)
    qc.rx(np.pi/2, 0)

    # Generate the process tomography circuits
    tomo_circuits = process_tomography_circuits(qc, [0])

    # Execute the tomography circuits on a simulator
    simulator = Aer.get_backend('qasm_simulator')
    job = execute(tomo_circuits, simulator, shots=10000)

    # Fit the measurement results to obtain the process matrix
    tomo_fitter = ProcessTomographyFitter(job.result(), tomo_circuits)
    choi_matrix = tomo_fitter.fit()

    # Print the estimated Choi matrix
    print(choi_matrix)

    # Plot the estimated Choi matrix as a cityscape
    plot_state_city(choi_matrix)


# Hadamard Test
def hadamard_test():
    # Create the circuit
    qc = QuantumCircuit(3, 1)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(0, 2)
    qc.h(0)
    qc.measure(0, 0)

    # Execute the circuit on a simulator
    simulator = Aer.get_backend('qasm_simulator')
    job = execute(qc, simulator, shots=1024)

    # Get the measurement result
    result = job.result()
    counts = result.get_counts(qc)
    print(counts)

    # Plot the Bloch vector of qubit 1
    backend = Aer.get_backend('statevector_simulator')
    job = execute(qc, backend)
    statevector = job.result().get_statevector(qc)
    plot_bloch_multivector(statevector)


# Main function to run the examples
def main():
    print("Quantum State Tomography:")
    quantum_state_tomography()

    print("\nQuantum Process Tomography:")
    quantum_process_tomography()