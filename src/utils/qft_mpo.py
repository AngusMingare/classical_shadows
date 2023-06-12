import numpy as np

def qft_mpo_3_qubits():

    ##############
    # Definitions
    ##############

    Id_gate = np.identity(2) 
    zero_gate = np.zeros((2,2))
    projector_00 = np.array([[1,0],[0,0]])
    projector_11 = np.array([[0,0],[0,1]])
    copy_gate = np.array([projector_00, projector_11])
    H_gate = np.array([[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),-1/np.sqrt(2)]])

    def phase_gate(theta):
        return np.array([[1,0],[0,np.exp(complex(0,1)*theta)]])

    # index order up, right, left
    def phase_tensor_rank_3(theta):
        return np.array([Id_gate, phase_gate(theta)])

    # index order up, down, right, left
    def phase_tensor_rank_4(theta):
        return np.array([[Id_gate, zero_gate], [zero_gate, phase_gate(theta)]])

    ####################
    # Zip-up Algorithm
    ####################

    # First contract all H gates into the neighbouring tensor
    contraction1 = np.einsum('dcb,ba->dca', copy_gate.copy(), H_gate.copy())
    contraction2 = np.einsum('gf,dife->dige', H_gate.copy(), phase_tensor_rank_4(np.pi/2))
    contraction3 = np.einsum('nm,jml->jnl', H_gate.copy(), phase_tensor_rank_3(np.pi/2))

    # Now start the zip-up process
    contraction4 = np.einsum('jnl,ilk->ijnk', contraction3, phase_tensor_rank_3(np.pi/4))
    contraction4 = np.reshape(contraction4, (4,4))
    u1, s1, vh1 = np.linalg.svd(contraction4, full_matrices=False)
    u1s1 = u1[:, :4] * s1
    u1s1 = np.reshape(u1s1, (2,2,4))
    vh1 = np.reshape(vh1, (4,2,2))

    contraction5 = np.einsum('jhg,dige->dijhe', copy_gate.copy(), contraction2)
    contraction6 = np.einsum('dijhe,ijo->dohe', contraction5, u1s1)

    # Move the orthogonality centre back to the bottom
    contraction6 = np.moveaxis(contraction6, [1], [0])
    contraction6 = np.reshape(contraction6, (4,8))
    u2, s2, vh2 = np.linalg.svd(contraction6, full_matrices=False)
    u2s2 = u2[:, :4] * s2
    u2s2 = u2s2.T
    vh2 = np.reshape(vh2, (4,2,2,2))
    vh2 = np.moveaxis(vh2, [0], [1])
    contraction7 = np.einsum('qo,onk->qnk', u2s2, vh1)

    # The final matrices in the MPO
    final1 = contraction1
    final2 = vh2 
    final3 = contraction7

    return final1, final2, final3

def qft_mpo_n_qubits_no_truncation(n):

    ##############
    # Definitions
    ##############

    Id_gate = np.identity(2) 
    zero_gate = np.zeros((2,2))
    projector_00 = np.array([[1,0],[0,0]])
    projector_11 = np.array([[0,0],[0,1]])
    copy_gate = np.array([projector_00, projector_11])
    H_gate = np.array([[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),-1/np.sqrt(2)]])

    def phase_gate(theta):
        return np.array([[1,0],[0,np.exp(complex(0,1)*theta)]])

    # index order up, right, left
    def phase_tensor_rank_3(theta):
        return np.array([Id_gate, phase_gate(theta)])

    # index order up, down, right, left
    def phase_tensor_rank_4(theta):
        return np.array([[Id_gate, zero_gate], [zero_gate, phase_gate(theta)]])
    
    if n == 3:
        return qft_mpo_3_qubits()
    
    else:
        # Contract the H gate
        tensor_1 = np.einsum('dcb,ba->dca', copy_gate.copy(), H_gate.copy())

        ##########
        # Zip up
        ##########
        qft_mpo_n_minus_1_qubits = qft_mpo_n_qubits_no_truncation(n-1)
        output_tensors = []

        # First Contraction
        final_n_minus_1 = qft_mpo_n_minus_1_qubits[-1]
        tensor_n = np.einsum('mkj,lji->lmki', final_n_minus_1, phase_tensor_rank_3(np.pi/(2**(n-1))))

        # First SVD
        tensor_n_rows = tensor_n.shape[0] * tensor_n.shape[1]
        tensor_n = np.reshape(tensor_n, (tensor_n_rows, 4))
        un, sn, tensor_n = np.linalg.svd(tensor_n, full_matrices=False)
        num_singular_vals = sn.shape[0]
        unsn = un[:, :num_singular_vals] * sn
        unsn = np.reshape(unsn, (2,final_n_minus_1.shape[0],num_singular_vals))
        tensor_n = np.reshape(tensor_n, (num_singular_vals,2,2))

        output_tensors.append(tensor_n)

        feed_in_matrix = unsn
        for val in range(2,n-1):
            # Contraction
            final_n_minus_val = qft_mpo_n_minus_1_qubits[-val]
            tensor_n_minus_val = np.einsum('vmkj,ulji->uvlmki', final_n_minus_val, phase_tensor_rank_4(np.pi/(2**(n-val))))
            tensor_n_minus_val = np.einsum('uvlmki,lms->uvski', tensor_n_minus_val, feed_in_matrix)

            # SVD
            tensor_n_minus_val = np.reshape(tensor_n_minus_val, (tensor_n_minus_val.shape[0]*tensor_n_minus_val.shape[1], tensor_n_minus_val.shape[2]*tensor_n_minus_val.shape[3]*tensor_n_minus_val.shape[4]))
            u_n_minus_val, s_n_minus_val, tensor_n_minus_val = np.linalg.svd(tensor_n_minus_val, full_matrices=False)
            num_singular_vals = s_n_minus_val.shape[0]
            u_n_minus_val_s_n_minus_val = u_n_minus_val[:, :num_singular_vals] * s_n_minus_val 
            u_n_minus_val_s_n_minus_val = np.reshape(u_n_minus_val_s_n_minus_val, (2, final_n_minus_val.shape[0], num_singular_vals))
            tensor_n_minus_val = np.reshape(tensor_n_minus_val, (num_singular_vals, feed_in_matrix.shape[2], 2, 2))

            output_tensors.append(tensor_n_minus_val)
            feed_in_matrix = u_n_minus_val_s_n_minus_val

        # Final Contraction
        final_1 = qft_mpo_n_minus_1_qubits[0]
        tensor_2 = np.einsum('mkj,dlji->dlmki', final_1, phase_tensor_rank_4(np.pi/2))
        tensor_2 = np.einsum('dlmki,lmn->dnki', tensor_2, feed_in_matrix)
        output_tensors.append(tensor_2)

        # Reverse the list order
        output_tensors = output_tensors[::-1]

        # Return the orthogonality centre to the bottom
        final_tensors = [tensor_1]

        tensor_2 = np.moveaxis(tensor_2, [1], [0])
        tensor_2 = np.reshape(tensor_2, (tensor_2.shape[0], tensor_2.shape[1]*tensor_2.shape[2]*tensor_2.shape[3]))
        u2, s2, tensor_2 = np.linalg.svd(tensor_2, full_matrices=False)
        num_singular_vals = s2.shape[0]
        u2s2 = u2[:, :num_singular_vals] * s2
        u2s2 = u2s2.T
        tensor_2 = np.reshape(tensor_2, (num_singular_vals, 2,2,2))
        tensor_2 = np.moveaxis(tensor_2, [0], [1])
        final_tensors.append(tensor_2)

        feed_in_matrix = u2s2
        for val in range(3,n):
            tensor_val = output_tensors[val-2]
            tensor_val = np.einsum('adcb,ea->edcb', tensor_val, feed_in_matrix)
            tensor_val = np.moveaxis(tensor_val, [1], [0])
            tensor_val = np.reshape(tensor_val, (tensor_val.shape[0], tensor_val.shape[1]*tensor_val.shape[2]*tensor_val.shape[3]))
            uval, sval, tensor_val = np.linalg.svd(tensor_val, full_matrices=False)
            num_singular_vals = sval.shape[0]
            uvalsval = uval[:, :num_singular_vals] * sval
            uvalsval = uvalsval.T
            tensor_val = np.reshape(tensor_val, (num_singular_vals,feed_in_matrix.shape[0],2,2))
            tensor_val = np.moveaxis(tensor_val, [0], [1])

            final_tensors.append(tensor_val)
            feed_in_matrix = uvalsval
        
        tensor_n = output_tensors[-1]
        tensor_n = np.einsum('ab,bdc->adc', feed_in_matrix, tensor_n)
        final_tensors.append(tensor_n)

        return final_tensors
    
def qft_mpo_n_qubits(n, max_bond_dimension):

    ##############
    # Definitions
    ##############

    Id_gate = np.identity(2) 
    zero_gate = np.zeros((2,2))
    projector_00 = np.array([[1,0],[0,0]])
    projector_11 = np.array([[0,0],[0,1]])
    copy_gate = np.array([projector_00, projector_11])
    H_gate = np.array([[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),-1/np.sqrt(2)]])

    def phase_gate(theta):
        return np.array([[1,0],[0,np.exp(complex(0,1)*theta)]])

    # index order up, right, left
    def phase_tensor_rank_3(theta):
        return np.array([Id_gate, phase_gate(theta)])

    # index order up, down, right, left
    def phase_tensor_rank_4(theta):
        return np.array([[Id_gate, zero_gate], [zero_gate, phase_gate(theta)]])
    
    if n == 3:
        return qft_mpo_3_qubits()
    
    else:
        # Contract the H gate
        tensor_1 = np.einsum('dcb,ba->dca', copy_gate.copy(), H_gate.copy())

        ##########
        # Zip up
        ##########
        qft_mpo_n_minus_1_qubits = qft_mpo_n_qubits_no_truncation(n-1)
        output_tensors = []

        # First Contraction
        final_n_minus_1 = qft_mpo_n_minus_1_qubits[-1]
        tensor_n = np.einsum('mkj,lji->lmki', final_n_minus_1, phase_tensor_rank_3(np.pi/(2**(n-1))))

        # First SVD
        tensor_n_rows = tensor_n.shape[0] * tensor_n.shape[1]
        tensor_n = np.reshape(tensor_n, (tensor_n_rows, 4))
        un, sn, tensor_n = np.linalg.svd(tensor_n, full_matrices=False)
        num_singular_vals = sn.shape[0]
        bond_dim = max_bond_dimension if max_bond_dimension < num_singular_vals else num_singular_vals
        unsn = un[:, :bond_dim] * sn[:bond_dim]
        unsn = np.reshape(unsn, (2,final_n_minus_1.shape[0],bond_dim))
        tensor_n = tensor_n[:bond_dim, :]
        tensor_n = np.reshape(tensor_n, (bond_dim,2,2))

        output_tensors.append(tensor_n)

        feed_in_matrix = unsn
        for val in range(2,n-1):
            # Contraction
            final_n_minus_val = qft_mpo_n_minus_1_qubits[-val]
            tensor_n_minus_val = np.einsum('vmkj,ulji->uvlmki', final_n_minus_val, phase_tensor_rank_4(np.pi/(2**(n-val))))
            tensor_n_minus_val = np.einsum('uvlmki,lms->uvski', tensor_n_minus_val, feed_in_matrix)

            # SVD
            tensor_n_minus_val = np.reshape(tensor_n_minus_val, (tensor_n_minus_val.shape[0]*tensor_n_minus_val.shape[1], tensor_n_minus_val.shape[2]*tensor_n_minus_val.shape[3]*tensor_n_minus_val.shape[4]))
            u_n_minus_val, s_n_minus_val, tensor_n_minus_val = np.linalg.svd(tensor_n_minus_val, full_matrices=False)
            num_singular_vals = s_n_minus_val.shape[0]
            bond_dim = max_bond_dimension if max_bond_dimension < num_singular_vals else num_singular_vals
            u_n_minus_val_s_n_minus_val = u_n_minus_val[:, :bond_dim] * s_n_minus_val[:bond_dim]
            u_n_minus_val_s_n_minus_val = np.reshape(u_n_minus_val_s_n_minus_val, (2, final_n_minus_val.shape[0], bond_dim))
            tensor_n_minus_val = tensor_n_minus_val[:bond_dim, :]
            tensor_n_minus_val = np.reshape(tensor_n_minus_val, (bond_dim, feed_in_matrix.shape[2], 2, 2))

            output_tensors.append(tensor_n_minus_val)
            feed_in_matrix = u_n_minus_val_s_n_minus_val

        # Final Contraction
        final_1 = qft_mpo_n_minus_1_qubits[0]
        tensor_2 = np.einsum('mkj,dlji->dlmki', final_1, phase_tensor_rank_4(np.pi/2))
        tensor_2 = np.einsum('dlmki,lmn->dnki', tensor_2, feed_in_matrix)
        output_tensors.append(tensor_2)

        # Reverse the list order
        output_tensors = output_tensors[::-1]

        # Return the orthogonality centre to the bottom
        final_tensors = [tensor_1]

        tensor_2 = np.moveaxis(tensor_2, [1], [0])
        tensor_2 = np.reshape(tensor_2, (tensor_2.shape[0], tensor_2.shape[1]*tensor_2.shape[2]*tensor_2.shape[3]))
        u2, s2, tensor_2 = np.linalg.svd(tensor_2, full_matrices=False)
        num_singular_vals = s2.shape[0]
        bond_dim = max_bond_dimension if max_bond_dimension < num_singular_vals else num_singular_vals
        u2s2 = u2[:, :bond_dim] * s2[:bond_dim]
        u2s2 = u2s2.T
        tensor_2 = tensor_2[:bond_dim, :]
        tensor_2 = np.reshape(tensor_2, (bond_dim, 2,2,2))
        tensor_2 = np.moveaxis(tensor_2, [0], [1])
        final_tensors.append(tensor_2)

        feed_in_matrix = u2s2
        for val in range(3,n):
            tensor_val = output_tensors[val-2]
            tensor_val = np.einsum('adcb,ea->edcb', tensor_val, feed_in_matrix)
            tensor_val = np.moveaxis(tensor_val, [1], [0])
            tensor_val = np.reshape(tensor_val, (tensor_val.shape[0], tensor_val.shape[1]*tensor_val.shape[2]*tensor_val.shape[3]))
            uval, sval, tensor_val = np.linalg.svd(tensor_val, full_matrices=False)
            num_singular_vals = sval.shape[0]
            bond_dim = max_bond_dimension if max_bond_dimension < num_singular_vals else num_singular_vals
            uvalsval = uval[:, :bond_dim] * sval[:bond_dim]
            uvalsval = uvalsval.T
            tensor_val = tensor_val[:bond_dim, :]
            tensor_val = np.reshape(tensor_val, (bond_dim,feed_in_matrix.shape[0],2,2))
            tensor_val = np.moveaxis(tensor_val, [0], [1])

            final_tensors.append(tensor_val)
            feed_in_matrix = uvalsval
        
        tensor_n = output_tensors[-1]
        tensor_n = np.einsum('ab,bdc->adc', feed_in_matrix, tensor_n)
        final_tensors.append(tensor_n)

        return final_tensors