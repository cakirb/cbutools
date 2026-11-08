#!/usr/bin/env ipython
"""Tests for the LARRYv2 recipe."""

import os.path
from pathlib import Path

import pandas as pd

import cbutools as cbu
from cbutools.recipes import LARRYV2_LIBRARIES

HERE = Path(__file__).parent

# The three motifs in the fixtures. Barcodes come back as the 6 nt library
# identifier followed by the 28 nt motif.
M1 = "AAAATGCCCCCATTTTGTGGGGAGACAC"
M2 = "TTAATGGGCCCAAATTGTCCGGAGTGTG"
M3 = "CACATGGAGACATCTCGTATATAGCTCT"
GFP = LARRYV2_LIBRARIES["egfp1"].identifier

CBC1, CBC2, CBC3 = "AACCTGAAGTATTCCG", "TTTCTGAAGTATTCCG", "GGGGTGAAGTATTCCG"
UMI1, UMI2, UMI3 = "CGATCAGATGAA", "CGATCAGATGGT", "CGATCAGATCCC"


def files(name):
    return {
        "r1": str(HERE / f"data/mock_larryv2_{name}_r1.fq"),
        "r2": str(HERE / f"data/mock_larryv2_{name}_r2.fq"),
    }


def test_cutadapt():
    """Checking if cutadapt extracts the 34 nt prefix+motif correctly

    The fixture's three good reads carry the 2-4 nt primer stagger, which the
    non-anchored 5' adapter has to absorb.
    """
    bcu = cbu.get_barcodes(files("cutadapt"), type="larryv2").sort_index()

    ind = pd.MultiIndex.from_tuples(
        [
            (CBC1, GFP + M1, UMI1),
            (CBC2, GFP + M2, UMI2),
            (CBC3, GFP + M3, UMI3),
        ]
    )
    ref = pd.Series([1, 1, 1], index=ind).sort_index()

    assert bcu.equals(ref)
    assert all(len(b) == 34 for b in bcu.index.get_level_values("Barcode"))


def test_wrong_spacers_are_rejected():
    """The fixed dinucleotides are a real gate, not decoration.

    The fixture's lv4_badspacers read has the right length and the right flanks
    but AA in place of all four spacers, so the motif stage must drop it.
    """
    bcu = cbu.get_barcodes(files("cutadapt"), type="larryv2")
    assert len(bcu) == 3
    for barcode in bcu.index.get_level_values("Barcode"):
        assert barcode[10:12] == "TG"
        assert barcode[16:18] == "CA"
        assert barcode[22:24] == "GT"
        assert barcode[28:30] == "AG"


def test_count():
    """Checking if reads are counted correctly"""
    bcu = cbu.get_barcodes(files("counting"), type="larryv2").sort_index()

    ind = pd.MultiIndex.from_tuples(
        [
            (CBC1, GFP + M1, UMI1),
            (CBC1, GFP + M1, UMI2),
            (CBC2, GFP + M2, UMI3),
            (CBC3, GFP + M3, UMI1),
        ]
    )
    ref = pd.Series([2, 1, 4, 5], index=ind).sort_index()

    assert bcu.equals(ref)


def test_count_CBC():
    """Checking the valid_CBC argument"""
    bcu = cbu.get_barcodes(
        files("counting"), type="larryv2", valid_CBC=[CBC3]
    ).sort_index()

    ind = pd.MultiIndex.from_tuples([(CBC3, GFP + M3, UMI1)])
    ref = pd.Series([5], index=ind).sort_index()

    assert bcu.equals(ref)


def test_every_library_is_kept_and_distinguishable():
    """One read per Addgene library in, five distinct barcodes out.

    All five reads share the same 28 nt motif, so the identifier is the only
    thing telling them apart. That is the point of keeping it: the same random
    bases under two libraries are two clones, and the library is recoverable
    from barcode[:6] with no second pass.
    """
    identifiers = {library.identifier for library in LARRYV2_LIBRARIES.values()}
    bcu = cbu.get_barcodes(files("libraries"), type="larryv2")

    barcodes = sorted(bcu.index.get_level_values("Barcode"))
    assert barcodes == sorted(identifier + M1 for identifier in identifiers)

    by_library = {b[:6]: b[6:] for b in barcodes}
    assert set(by_library) == identifiers
    assert set(by_library.values()) == {M1}


def test_process_larryv2_matches_dispatch():
    direct = cbu.larryv2.process_larryv2(files("libraries")).sort_index()
    dispatched = cbu.get_barcodes(files("libraries"), type="larryv2").sort_index()
    assert direct.equals(dispatched)


def test_debug():
    """Checking if debug=True works and files are available for inspection"""
    bcu, tempDIR = cbu.larryv2.process_larryv2(files("counting"), debug=True)

    assert isinstance(bcu, cbu.CBUclasses.CBUSeries)
    for number in (1, 2):
        assert os.path.isfile(f"{tempDIR.name}/bar{number}.fq")
        assert os.path.isfile(f"{tempDIR.name}/cbcumi{number}.fq")
    assert not os.path.isfile(tempDIR.name + "/bar3.fq")


def test_explicit_intermediate_dir(tmp_path):
    """Checking caller-controlled intermediate directories used by workflows."""
    intermediate_dir = tmp_path / ".intermediate_cutadapt"

    bcu, returned_dir = cbu.larryv2.process_larryv2(
        files("counting"),
        debug=True,
        cores=1,
        intermediate_dir=intermediate_dir,
    )

    assert isinstance(bcu, cbu.CBUclasses.CBUSeries)
    assert returned_dir == intermediate_dir.resolve()
    assert (returned_dir / "bar2.fq").is_file()
