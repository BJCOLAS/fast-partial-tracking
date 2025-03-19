import numpy as np


def time_vector_computation(n,Q) : 
    """
    Compute the time vector 

    Input : n is the discrete time instant
            Q is the polynomial order 

    Output : the time vector

    """
    return n ** np.arange(Q+1)
    

def derivative_time_vector_computation(n,Q):
    """
    Compute the derivative of the time vector 

    Input : n is the discrete time instant
           Q is the polynomial order 

    Output : the derivative time vector

    """
    return n ** np.arange(0,Q) * np.arange(1, Q+1)

def FreqAmpPha(fs, alpha, p, pPrime) :
    """
    Compute the log-amplitude, the frequency and the phase from alpha

    Input : 


    Output :

    """
    # Eq. (7)

    freq = np.array(np.expand_dims(fs/(2*np.pi) * np.imag(pPrime @ alpha[1:, :]),axis=1))

    # Eq. (6)

    Pha = np.array(np.expand_dims(np.imag(p @ alpha),axis=1))

    # Eq. (5)

    amp  = np.array(np.expand_dims(np.real(p @ alpha),axis=1))

    return np.array([freq, amp, Pha]).T[:]


def match_sets(A, B):
        """
        Matches the elements of A and B.
        
        Parameters:
        A (list or array): First set of unique elements.
        B (list or array): Second set of unique elements.
        
        Returns:
        tuple: (match, nomatch_A, nomatch_B)
            match (numpy array): Pairs of matching indices.
            nomatch_A (numpy array): Indices of elements in A that did not match.
            nomatch_B (numpy array): Indices of elements in B that did not match.
        """

        nA = len(A)
        nB = len(B)
        
        # Ensure elements in A and B are unique
        if len(set(A)) != nA:
            print("A must be unique")
            return None, None, None
        elif len(set(B)) != nB:
            print("B must be unique")
            return None, None, None
        
        validA = np.ones(nA, dtype=bool)  # Boolean mask for valid elements in A
        validB = np.ones(nB, dtype=bool)  # Boolean mask for valid elements in B
        
        match = np.zeros((min(nA, nB), 2), dtype=int)  # Array to store matched indices
        counter = 0  # Counter for matched pairs
        
        for a in range(nA):
            if validA[a]:
                # Find indices in B that match A[a] and are still valid
                indices = np.where(B[validB] == A[a])[0]
                
                if indices.size > 0:
                    counter += 1
                    bs = np.where(validB)[0]  # Get valid indices in B
                    b = bs[indices[0]]  # Select the first match
                    
                    match[counter - 1, :] = [a, b]
                    validA[a] = False
                    validB[b] = False
        
        match = match[:counter, :]  # Trim to the actual number of matches
        nomatch_A = np.where(validA)[0]  # Indices of unmatched elements in A
        nomatch_B = np.where(validB)[0]  # Indices of unmatched elements in B
        
        return match, nomatch_A, nomatch_B