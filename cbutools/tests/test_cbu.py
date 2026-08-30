import cbutools as cbu
import numpy as np
from pathlib import Path
import pandas as pd
import pytest

HERE = Path(__file__).parent


@pytest.fixture
def cbu_example():
    """A small example of CBUSeries save"""
    x = cbu.load_barcodes(HERE / "data/cbu_raw.csv")
    return x


@pytest.fixture
def cb_example():
    """A small example of CBSeries save"""
    x = cbu.load_barcodes(HERE / "data/cbu_umicount.csv")
    return x


def test_read_filter():
    pre = cbu.load_barcodes(HERE / "data/cbu_raw.csv")
    ref = cbu.load_barcodes(HERE / "data/cbu_readfiltered.csv")

    totest = pre.filter_by_reads(min_counts=4)
    assert ref.equals(totest)


def test_hamming_filter():
    pre = cbu.load_barcodes(HERE / "data/cbu_readfiltered.csv").sort_index()
    ref = cbu.load_barcodes(HERE / "data/cbu_hamfiltered.csv").sort_index()

    totest = pre.filter_by_hamming(min_distance=3, plot=False).sort_index()

    assert ref.equals(totest)
    assert len(totest) - len(pre) == -2


def test_count_umis():
    pre = cbu.load_barcodes(HERE / "data/cbu_hamfiltered.csv").sort_index()
    ref = cbu.load_barcodes(HERE / "data/cbu_umicount.csv").sort_index()

    totest = pre.count_UMI()
    print(totest)
    print(ref)

    assert ref.equals(totest)


def test_UMI_filter():

    pre = cbu.load_barcodes(HERE / "data/cbu_umicount.csv").sort_index()
    ref = cbu.load_barcodes(HERE / "data/cbu_umifiltered.csv").sort_index()

    totest = pre.filter_by_UMI(min_counts=2)
    assert ref.equals(totest)


def test_assign_barcodes():

    pre = cbu.load_barcodes(HERE / "data/cbu_umifiltered.csv").sort_index()
    ref = pd.read_csv(HERE / "data/cbu_assigned.csv")
    totest = pre.assign_barcodes()

    totest = totest.loc[ref.CBC, :]
    assert ref.CBC.tolist() == totest.index.tolist()
    assert ref.assigned.tolist() == totest.Barcode.tolist()


def test_assign_barcodes_dispr():
    """This is the same test as test_assign_barcodes but checks also the dispr_filter argument"""
    pre = cbu.load_barcodes(HERE / "data/cbu_umifiltered.csv").sort_index()
    ref = pd.read_csv(HERE / "data/cbu_assigned_dispr.csv")
    totest = pre.assign_barcodes(dispr_filter=0.51)

    totest = totest.loc[ref.CBC, :]
    assert ref.CBC.tolist() == totest.index.tolist()
    assert ref.assigned.tolist() == totest.Barcode.tolist()


def test_check_CBU_index(cbu_example):
    cbu_example.filter_by_reads()

    with pytest.raises(Exception):
        cbu_example.index.names = ["A", "B", "C"]
        cbu_example.filter_by_reads()


def test_check_CB_index(cb_example):
    cb_example.filter_by_UMI()

    with pytest.raises(Exception):
        cbu_example.index.names = ["A", "B"]
        cbu_example.filter_by_UMI()
