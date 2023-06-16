def encodeJordanWigner():
    """Jordan-Wigner transformation
    
    The Jordan-Wigner transformation maps a fermionic Hamiltonian
    to a spin (qubit) Hamiltonian for simulation on quantum computers.
    
    Parameters
    -----------
    hamiltonian : dict
        A dictionary of pairs {"indices" : weight} encoding the second quantised Hamiltonian. 
        The keys are assumed to be strings of indices of even length, the first n/2 corresponding to
        creation operators and the last n/2 corresponding to annihilation operators. 

        E.g., H = 2 a_0^dagger a_1 - 4 a_0^dagger a_1^dagger a_2 a_3 is written as
        {"01" : 2, "0123" : -4}
    
    Returns
    --------
    encoded_hamiltonian : dict
        Pairs {"pauli_string" : weight} for the encoded Hamiltonian 
    """

    return

def encodeBravyiKitaev():
    return 

def encodeParity():
    return

def qubitTapering():
    return 

def contextualSubspace():
    return