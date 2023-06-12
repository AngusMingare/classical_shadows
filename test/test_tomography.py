from ..src.tomography.qst import *
from ..src.tomography.qpt import *
from qiskit import QuantumCircuit
import numpy as np

class TestQST:
    def test_simple(self):
        prep_circ = QuantumCircuit(1)
        prep_circ.h(0)
        state_dm = quantumStateTomography(prep_circ, 0)

        assert np.round(state_dm[0][0] - 0.5, 1) == 0
        assert np.round(state_dm[0][1] - 0.5, 1) == 0
        assert np.round(state_dm[1][0] - 0.5, 1) == 0
        assert np.round(state_dm[1][1] - 0.5, 1) == 0

        return
    
    def test_entangled(self):
        prep_circ = QuantumCircuit(2)
        prep_circ.h(0)
        prep_circ.cnot(0,1)
        state_dm_0 = quantumStateTomography(prep_circ, 0)
        state_dm_1 = quantumStateTomography(prep_circ, 1)

        assert np.round(state_dm_0[0][0] - 0.5, 1) == 0
        assert np.round(state_dm_0[0][1], 1) == 0
        assert np.round(state_dm_0[1][0], 1) == 0
        assert np.round(state_dm_0[1][1] - 0.5, 1) == 0

        assert np.round(state_dm_1[0][0] - 0.5, 1) == 0
        assert np.round(state_dm_1[0][1], 1) == 0
        assert np.round(state_dm_1[1][0], 1) == 0
        assert np.round(state_dm_1[1][1] - 0.5, 1) == 0

        return
    
class TestQPT:
    def test_simple(self):
        circ = QuantumCircuit(1)
        circ.h(0)
        chi = quantumProcessTomography(circ, 0)
        choi = chiToChoi(chi)
        ev = choiToEvolution(choi)

        assert all([np.round(np.abs(np.real(ev[m,n])) - 0.5, 1) == 0 for m in range(4) for n in range(4)])

        return