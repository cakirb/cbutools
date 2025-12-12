import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def sequences_to_numeric(seqs: np.ndarray) -> np.ndarray:
    """Convert string sequences to numeric matrix for vectorized operations.
    
    Parameters
    ----------
    seqs : np.ndarray
        1D array of equal-length string sequences
    
    Returns
    -------
    np.ndarray
        2D array of shape (n_sequences, seq_length) with dtype=np.uint8
    """
    if len(seqs) == 0:
        return np.array([], dtype=np.uint8).reshape(0, 0)
    
    # Convert strings to 2D char array, then to numeric
    char_array = np.array([list(s) for s in seqs], dtype='U1')
    return char_array.view(np.uint32).astype(np.uint8)


def compute_hamming_matrix(seqs: np.ndarray) -> np.ndarray:
    """Compute pairwise Hamming distance matrix using vectorized operations.
    
    Parameters
    ----------
    seqs : np.ndarray
        1D array of equal-length string sequences
    
    Returns
    -------
    np.ndarray
        Symmetric n×n matrix of pairwise Hamming distances (dtype=np.int32)
    """
    if len(seqs) == 0:
        return np.array([], dtype=np.int32).reshape(0, 0)
    
    if len(seqs) == 1:
        return np.zeros((1, 1), dtype=np.int32)
    
    # Convert to numeric matrix (n_seqs × seq_length)
    numeric = sequences_to_numeric(seqs)
    
    # Compute all pairwise distances using broadcasting
    # Shape: (n, 1, seq_len) vs (1, n, seq_len) → (n, n, seq_len)
    # Sum mismatches along last axis → (n, n)
    hdist = (numeric[:, np.newaxis, :] != numeric[np.newaxis, :, :]).sum(axis=2)
    
    return hdist.astype(np.int32)


def get_sequence_length(seqs: np.ndarray) -> int | None:
    """Get the length of sequences in the array.
    
    Parameters
    ----------
    seqs : np.ndarray
        1D array of string sequences
    
    Returns
    -------
    int or None
        Length of sequences, or None if empty
    """
    if len(seqs) == 0:
        return None
    return len(seqs[0])


def hamming_filter(counts, min_distance=3):
    """Filters out barcodes based on Hamming distance.
    
    For each pair of barcodes with a Hamming distance below the threshold,
    the barcode with the lower count is rejected. Ties are recorded.
    
    Parameters
    ----------
    counts : pd.Series
        Pandas Series with sequences as the index and counts as values
    min_distance : int, default 3
        Minimum Hamming distance threshold
    
    Returns
    -------
    tokeep : pd.Index
        Index of sequences that passed the filter
    toreject : list
        List of sequences to be rejected
    ties : dict
        Dictionary mapping a sequence to tied sequences
    """
    print(f"Filtering using Hamming distance threshold of {min_distance}...")
    seqs = counts.index.values.astype(str)
    n = len(seqs)
    
    if n == 0:
        return counts.index, [], {}
    
    # Compute pairwise Hamming distances
    hdist = compute_hamming_matrix(seqs)
    
    # Plot histogram
    if n > 1:
        triu_indices = np.triu_indices_from(hdist, k=1)
        plt.hist(hdist[triu_indices], bins=50)
        plt.title("Pairwise Hamming Distance")
        plt.xlabel("Hamming Distance")
        plt.ylabel("Frequency")
        plt.show()
    
    toreject = set()
    ties = {}
    counts_arr = counts.values
    
    for i in range(n):
        seq = seqs[i]
        if seq in toreject:
            continue
        
        xcount = counts_arr[i]
        row = hdist[i, :]
        
        # Find neighbors below threshold (excluding self)
        mask = (row < min_distance) & (row > 0)
        if not mask.any():
            continue
        
        belowtr_indices = np.where(mask)[0]
        belowtr_counts = counts_arr[belowtr_indices]
        
        if np.any(xcount < belowtr_counts):
            toreject.add(seq)
        elif xcount == belowtr_counts.max():
            tied_mask = belowtr_counts == belowtr_counts.max()
            ties[seq] = list(seqs[belowtr_indices[tied_mask]])
    
    tokeep = counts.index[~counts.index.isin(toreject)]
    
    print(f"Passed: {len(tokeep)}\nRejected: {len(toreject)}\nTies: {len(ties)}")
    
    return tokeep, list(toreject), ties