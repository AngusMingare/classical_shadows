class SuperconductingDevice:

    def __init__(self, 
                 num_qubits : int = 16,
                 connectivity : str = "square",
                 basis_gates : list[str] = ["rz", "sx", "x", "id", "ecr"],
                 t1 : float = 100e-6,
                 t2 : float = 100e-6,
                 one_qubit_gate_time : float = 100e-9,
                 two_qubit_gate_time : float = 500e-9,
                 one_qubit_gate_fidelity : float = 1e-4,
                 two_qubit_gate_fidelity : float = 1e-2,
                 initialisation_time : float = 250e-6,
                 readout_time : float = 5e-6,
                 spam_fidelity : float = 1e-2) -> None:

        self.t1 = t1
        self.t2 = t2
        self.one_qubit_gate_time = one_qubit_gate_time
        self.two_qubit_gate_time = two_qubit_gate_time
        self.one_qubit_gate_fidelity = one_qubit_gate_fidelity
        self.two_qubit_gate_fidelity = two_qubit_gate_fidelity
        self.initialisation_time = initialisation_time
        self.readout_time = readout_time
        self.spam_fidelity = spam_fidelity

        if connectivity == "square":
            self.connectivity = "square"
            self.connectivity_graph = []

class SemiconductingDevice:

    def __init__(self, 
                 num_qubits : int = 16,
                 connectivity : str = "square",
                 basis_gates : list[str] = [],
                 t1 : float = 1,
                 t2 : float = 1,
                 one_qubit_gate_time : float = 1,
                 two_qubit_gate_time : float = 1,
                 one_qubit_gate_fidelity : float = 1,
                 two_qubit_gate_fidelity : float = 1,
                 initialisation_time : float = 1,
                 readout_time : float = 1,
                 spam_fidelity : float = 1) -> None:

        self.t1 = t1
        self.t2 = t2
        self.one_qubit_gate_time = one_qubit_gate_time
        self.two_qubit_gate_time = two_qubit_gate_time
        self.one_qubit_gate_fidelity = one_qubit_gate_fidelity
        self.two_qubit_gate_fidelity = two_qubit_gate_fidelity
        self.initialisation_time = initialisation_time
        self.readout_time = readout_time
        self.spam_fidelity = spam_fidelity

class TrappedIonDevice:

    def __init__(self, 
                 num_qubits : int = 16,
                 connectivity : str = "all-to-all",
                 basis_gates : list[str] = ["rxy", "rz", "zz"],
                 t1 : float = 10,
                 t2 : float = 1,
                 one_qubit_gate_time : float = 1e-5,
                 two_qubit_gate_time : float = 1e-3,
                 one_qubit_gate_fidelity : float = 5e-5,
                 two_qubit_gate_fidelity : float = 2e-3,
                 initialisation_time : float = 3e-5,
                 readout_time : float = 2e-4,
                 spam_fidelity : float = 1e-3) -> None:

        self.t1 = t1
        self.t2 = t2
        self.one_qubit_gate_time = one_qubit_gate_time
        self.two_qubit_gate_time = two_qubit_gate_time
        self.one_qubit_gate_fidelity = one_qubit_gate_fidelity
        self.two_qubit_gate_fidelity = two_qubit_gate_fidelity
        self.initialisation_time = initialisation_time
        self.readout_time = readout_time
        self.spam_fidelity = spam_fidelity

class DiamondDevice:

    def __init__(self, 
                 num_qubits : int = 16,
                 connectivity : str = "square",
                 basis_gates : list[str] = ["rx", "ry", "cz"],
                 t1 : float = 1,
                 t2 : float = 1e-3,
                 one_qubit_gate_time : float = 1e-6,
                 two_qubit_gate_time : float = 1e-6,
                 one_qubit_gate_fidelity : float = 1e-2,
                 two_qubit_gate_fidelity : float = 1e-1,
                 initialisation_time : float = 1,
                 readout_time : float = 1,
                 spam_fidelity : float = 1e-1) -> None:

        self.t1 = t1
        self.t2 = t2
        self.one_qubit_gate_time = one_qubit_gate_time
        self.two_qubit_gate_time = two_qubit_gate_time
        self.one_qubit_gate_fidelity = one_qubit_gate_fidelity
        self.two_qubit_gate_fidelity = two_qubit_gate_fidelity
        self.initialisation_time = initialisation_time
        self.readout_time = readout_time
        self.spam_fidelity = spam_fidelity
