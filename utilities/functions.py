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

def bh_window(N, d=1):
    """
    Blackman-Harris Window
    
    Parameters:
    N : int
        Length of window
    d : int, optional
        Order of derivative (default is 1)
        
    Returns:
    win : ndarray
        The computed window
    n : ndarray
        The time indices
    """
    if N % 2:
        Nwinhf = (N - 1) // 2
        n = np.arange(-Nwinhf, Nwinhf + 1)
    else:
        Nwinhf = N // 2 - 1
        n = np.arange(-Nwinhf - 1, Nwinhf + 1)
    
    a = np.array([0.35875, 0.48829, 0.14128, 0.01168])
    in_vals = np.array([2, 4, 6]) * np.pi / N
    
    if d == 1:
        win = a[0] + a[1] * np.cos(in_vals[0] * n) + a[2] * np.cos(in_vals[1] * n) + a[3] * np.cos(in_vals[2] * n)
    elif d == 2:
        win = -a[1] * in_vals[0] * np.sin(in_vals[0] * n) \
              -a[2] * in_vals[1] * np.sin(in_vals[1] * n) \
              -a[3] * in_vals[2] * np.sin(in_vals[2] * n)
    elif d == 3:
        win = -a[1] * in_vals[0]**2 * np.cos(in_vals[0] * n) \
              -a[2] * in_vals[1]**2 * np.cos(in_vals[1] * n) \
              -a[3] * in_vals[2]**2 * np.cos(in_vals[2] * n)
    else:
        raise ValueError("Derivative order d must be 1, 2, or 3")
    
    return win

def hann_window(N, derivative=1):
    """
    Hann Window
    
    Parameters:
    N : int
        Length of window
    derivative : int, optional
        Order of derivative (default is 1)
        
    Returns:
    win : ndarray
        The computed window
    """
    if N % 2 == 1:
        Nhf = (N - 1) // 2
        n = np.arange(0, N) - Nhf
        in_val = 2 * np.pi / (N - 1)
        
        if derivative == 1:
            win = 0.5 + 0.5 * np.cos(in_val * n)
        elif derivative == 2:
            win = -0.5 * in_val * np.sin(in_val * n)
        elif derivative == 3:
            win = -0.5 * in_val**2 * np.cos(in_val * n)
        else:
            raise ValueError("Derivative order must be 1, 2, or 3")
    else:
        n = np.arange(0, N)
        in_val = 2 * np.pi / (N - 1)
        
        if derivative == 1:
            win = 0.5 - 0.5 * np.cos(in_val * n)
        elif derivative == 2:
            win = 0.5 * in_val * np.sin(in_val * n)
        elif derivative == 3:
            win = 0.5 * in_val**2 * np.cos(in_val * n)
        else:
            raise ValueError("Derivative order must be 1, 2, or 3")
    
    return win
