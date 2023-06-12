from .qpt import *
import numpy as np
from qiskit import QuantumCircuit, Aer, execute
from qiskit.extensions import *

def quantumAssistedPhaseEstimation(unitary, num_precision_bits):
    myGate = UnitaryGate(unitary)
    c_myGate = myGate.control(1)

    circ = QuantumCircuit(1+int(np.log2(unitary.shape[0])))
    circ.append(c_myGate)

    chi = quantumProcessTomography(circ, 0)
    choi = chiToChoi(chi)
    ev = choiToEvolution(choi)


    return 