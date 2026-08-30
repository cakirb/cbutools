#!/usr/bin/env ipython

import cbutools as cbu
import numpy as np
import pandas as pd
import pytest


def test_compute_hamming_matrix():
    seqs = np.array(["AAAA", "AAAT", "AATT", "TTTT"])
    observed = cbu.hamming.compute_hamming_matrix(seqs)
    expected = np.array(
        [
            [0, 1, 2, 4],
            [1, 0, 1, 3],
            [2, 1, 0, 2],
            [4, 3, 2, 0],
        ]
    )
    assert np.array_equal(observed, expected)


def test_compute_hamming_matrix_empty_and_single():
    assert cbu.hamming.compute_hamming_matrix(np.array([])).shape == (0, 0)
    assert np.array_equal(
        cbu.hamming.compute_hamming_matrix(np.array(["AAAA"])),
        np.zeros((1, 1), dtype=np.int32),
    )


def test_compute_hamming_matrix_rejects_unequal_lengths():
    with pytest.raises(ValueError, match="unequal length"):
        cbu.hamming.compute_hamming_matrix(np.array(["AAAA", "AAA"]))


def test_hamming():
    x = pd.Series(dict(AAAAAAA=20, AAAABBB=10, BAAAAAB=8, CCCCCCC=10, CCDDCCC=10))

    keep, reject, tie = cbu.hamming.hamming_filter(
        x, min_distance=4, plot=False
    )

    assert (keep == np.array(["AAAAAAA", "CCCCCCC", "CCDDCCC"])).all()
    assert (reject == np.array(["AAAABBB", "BAAAAAB"])).all()
    assert tie == dict(CCCCCCC=["CCDDCCC"], CCDDCCC=["CCCCCCC"])
