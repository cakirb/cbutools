import matplotlib.pyplot as plt
import numpy as np


def sequences_to_numeric(seqs):
    """Convert equal-length string sequences to a numeric matrix."""
    if len(seqs) == 0:
        return np.empty((0, 0), dtype=np.uint8)

    if len({len(seq) for seq in seqs}) != 1:
        raise ValueError("Undefined for sequences of unequal length")

    char_array = np.array([list(seq) for seq in seqs], dtype="U1")
    return char_array.view(np.uint32).astype(np.uint8)


def compute_hamming_matrix(seqs):
    """Compute a full matrix of pairwise Hamming distances."""
    if len(seqs) == 0:
        return np.empty((0, 0), dtype=np.int32)
    if len(seqs) == 1:
        return np.zeros((1, 1), dtype=np.int32)

    numeric = sequences_to_numeric(seqs)
    n, sequence_length = numeric.shape
    distances = np.empty((n, n), dtype=np.int32)

    # Bound the temporary broadcast array while retaining the public square
    # matrix returned by the original implementation.
    max_temporary_elements = 64 * 1024 * 1024
    rows_per_block = max(
        1,
        min(
            n,
            max_temporary_elements // max(1, n * sequence_length),
        ),
    )
    for start in range(0, n, rows_per_block):
        stop = min(start + rows_per_block, n)
        distances[start:stop] = np.sum(
            numeric[start:stop, np.newaxis, :] != numeric[np.newaxis, :, :],
            axis=2,
            dtype=np.int32,
        )

    return distances


def hamming_filter(counts, min_distance=3, plot=True):
    """Filter less-abundant sequences below a Hamming-distance threshold.

    Parameters
    ----------
    counts : pandas.Series
        Sequences in the index and counts in the values.
    min_distance : int
        Sequences with distance below this threshold are compared.
    plot : bool
        Whether to display a histogram of pairwise Hamming distances.

    Returns
    -------
    pandas.Index, list, dict
        Sequences to keep, sequences to reject, and tied sequences.
    """
    print(f"Filtering using Hamming distance threshold of {min_distance}...")
    seqs = counts.index.values.astype(str)
    n = len(seqs)

    if n == 0:
        return counts.index, [], {}

    hdist = compute_hamming_matrix(seqs)

    if plot and n > 1:
        upper_triangle = np.triu_indices_from(hdist, k=1)
        plt.hist(hdist[upper_triangle], bins=50)
        plt.title("Pairwise Hamming Distance")
        plt.xlabel("Hamming Distance")
        plt.ylabel("Frequency")
        plt.show()

    toreject = []
    ties = {}
    count_values = counts.values

    for i, seq in enumerate(seqs):
        below_threshold = (hdist[i] < min_distance) & (np.arange(n) != i)
        if not below_threshold.any():
            continue

        neighbor_indices = np.flatnonzero(below_threshold)
        neighbor_counts = count_values[neighbor_indices]
        seq_count = count_values[i]

        if np.any(seq_count < neighbor_counts):
            toreject.append(seq)
        elif seq_count == neighbor_counts.max():
            tied = neighbor_indices[neighbor_counts == neighbor_counts.max()]
            ties[seq] = list(seqs[tied])

    tokeep = counts.index[~counts.index.isin(toreject)]
    print(f"Passed: {len(tokeep)}\nRejected: {len(toreject)}\nTies: {len(ties)}")
    return tokeep, toreject, ties
