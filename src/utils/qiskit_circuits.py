from qiskit import QuantumCircuit
import numpy as np

def qft_circ(n, measure=True):
    circ = QuantumCircuit(n,n)
    for control_index in range(n):
        circ.h(control_index)
        for rot_power in range(1, n-control_index):
            target_index = control_index+rot_power
            circ.cp(np.pi/(2**(rot_power)), control_index, target_index)
        circ.barrier()
    if measure:
        circ.measure(range(n),range(n))
    return circ

def iqft_circ(n, measure=True):
    circ = QuantumCircuit(n,n)
    for control_index in list(range(n))[::-1]:
        for rot_power in list(range(1,n-control_index))[::-1]:
            target_index = control_index+rot_power
            circ.cp(-np.pi/(2**(rot_power)), control_index, target_index)
        circ.h(control_index)
        circ.barrier()
    if measure:
        circ.measure(range(n), range(n))
    return circ
