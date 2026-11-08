#!/usr/bin/env ipython
"""Tests for the recipe registry and the metadata it exposes."""

import json
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import pytest

import cbutools as cbu
from cbutools.extract import process_recipe
from cbutools.recipes import (
    LARRYV2_LIBRARIES,
    LARRYV2_MIN_PREFIX_DISTANCE,
    CutadaptStage,
    Library,
    Recipe,
    Thresholds,
    get_recipe,
    list_recipes,
    recipe_metadata,
    register_recipe,
)

HERE = Path(__file__).parent
REPO_ROOT = HERE.parent.parent


RECIPE_NAMES = ["larry", "larryv2"]


def test_registered_recipes():
    assert list_recipes() == RECIPE_NAMES
    # Nothing is pending at the moment, so the two listings agree.
    assert list_recipes(include_pending=True) == RECIPE_NAMES


def test_unknown_recipe_lists_alternatives():
    with pytest.raises(KeyError) as excinfo:
        get_recipe("nosuchthing")
    message = str(excinfo.value)
    assert "nosuchthing" in message
    assert "larry" in message


def test_larry_chemistry_unchanged():
    """Regression guard: the refactor must not have altered LARRY extraction.

    These are the exact adapters process_larry used before recipes existed.
    """
    larry = get_recipe("larry")
    assert larry.barcode_length == 40
    assert (larry.cbc_length, larry.umi_length, larry.cbcumi_length) == (16, 12, 28)
    assert [(s.adapter, s.action) for s in larry.stages] == [
        ("TTGCTAGGAGAGACCATATG...ATGTCTGGATCCGATATCGC", None),
        ("NNNNTGNNNNCANNNNACNNNNGANNNNGTNNNNAGNNNN;min_overlap=36", "retain"),
    ]


def test_split_cbcumi():
    larry = get_recipe("larry")
    cbc, umi = larry.split_cbcumi("AACCTGAAGTATTCCG" + "CGATCAGATGAA")
    assert cbc == "AACCTGAAGTATTCCG"
    assert umi == "CGATCAGATGAA"


def test_accepts_barcode():
    larry = get_recipe("larry")
    assert larry.accepts_barcode("A" * 40)
    assert not larry.accepts_barcode("A" * 39)


def test_metadata_is_json_serialisable():
    """nf-larry reads this; it has to survive a round trip through JSON."""
    everything = recipe_metadata()
    restored = json.loads(json.dumps(everything))
    assert set(restored) == set(list_recipes(include_pending=True))

    larry = restored["larry"]
    assert larry["collapse_metric"] == "hamming"
    assert larry["barcode_length"] == 40
    assert larry["defaults"] == {
        "min_reads": 8,
        "min_umis": 3,
        "collapse_distance": 3,
        "dispr_filter": 0.1,
    }
    assert larry["lenient"] == {
        "min_reads": 1,
        "min_umis": 2,
        "collapse_distance": 3,
        "dispr_filter": 0.1,
    }
    assert larry["runnable"] is True


def test_metadata_for_single_recipe():
    assert recipe_metadata("larryv2")["label"] == "LARRYv2"


def test_pending_recipe_is_described_but_not_runnable():
    """A recipe may be registered for its metadata before it can be run.

    Nothing in the package is pending now that LARRYv2 has its flanks, but the
    mechanism is how a half-known chemistry gets recorded without guessing.
    """
    pending = Recipe(
        name="notyet",
        label="Not Yet",
        description="",
        barcode_length=20,
        pending="The 3' flank is unknown.",
    )
    assert pending.runnable is False
    assert pending.to_dict()["pending"] == "The 3' flank is unknown."
    assert "notyet" not in list_recipes()

    with pytest.raises(NotImplementedError, match="not runnable yet"):
        process_recipe(pending, {"r1": "unused.fq", "r2": "unused.fq"})


def _hamming(a, b):
    assert len(a) == len(b)
    return sum(x != y for x, y in zip(a, b))


def test_library_identifiers_are_far_enough_apart_to_collapse_safely():
    """larryv2 must not be able to fuse barcodes from different libraries.

    It keeps the 6 nt library identifier inside a 34 bp barcode, and collapse
    runs on the whole string, so the identifier is treated as mutable sequence
    when it is not. A cross-library pair is at least LARRYV2_MIN_PREFIX_DISTANCE
    apart; collapsing is exclusive, so it is safe exactly while

        collapse_distance <= min pairwise identifier distance.

    This guards a bug that was real at collapse_distance=5: the closest
    identifiers (TGCTAA tsapphire2, TGATTG egfp2) are 3 apart, so two barcodes
    from different libraries whose random portions differed by one base would
    merge into one clone.

    It also guards the future: adding a library whose identifier sits closer to
    an existing one fails here rather than silently in someone's clone calls.
    """
    identifiers = [library.identifier for library in LARRYV2_LIBRARIES.values()]
    distances = [_hamming(a, b) for a, b in combinations(identifiers, 2)]
    assert min(distances) == LARRYV2_MIN_PREFIX_DISTANCE

    larryv2 = get_recipe("larryv2")
    for thresholds in (larryv2.defaults, larryv2.lenient):
        assert thresholds.collapse_distance <= min(distances)


def test_larryv2_recipe():
    """One recipe, prefix kept, cloneRanger's read/UMI cutoffs, safe collapse."""
    larryv2 = get_recipe("larryv2")
    assert larryv2.barcode_length == 34
    assert len(larryv2.stages) == 2
    assert larryv2.stages[1].adapter.startswith("NNNNNN")
    assert larryv2.stages[1].adapter.endswith(";e=0")
    # reads_cutoff: 5 and umi_cutoff: 3 from cloneRanger's larry_config.yaml.
    assert larryv2.defaults.min_reads == 5
    assert larryv2.defaults.min_umis == 3
    # Not cloneRanger's d<=4 (which would be 5): see the invariant test above.
    assert larryv2.defaults.collapse_distance == 3
    assert larryv2.lenient.collapse_distance == 3


@pytest.mark.parametrize("length", [None, 0])
def test_recipe_requires_a_fixed_barcode_length(length):
    with pytest.raises(ValueError, match="fixed barcode_length"):
        Recipe(
            name="bad",
            label="Bad",
            description="",
            barcode_length=length,
            stages=(CutadaptStage(adapter="A...B"),),
        )


def test_larryv2_metadata_lists_the_addgene_libraries():
    """nf-larry names libraries from this metadata instead of copying them."""
    libraries = recipe_metadata("larryv2")["libraries"]
    assert {key: value["identifier"] for key, value in libraries.items()} == {
        "tsapphire1": "GTTCCA",
        "tsapphire2": "TGCTAA",
        "egfp1": "TCCAGT",
        "egfp2": "TGATTG",
        "mscarlet": "AAGATT",
    }
    assert libraries["tsapphire1"] == {
        "identifier": "GTTCCA",
        "label": "LARRYv2-T-Sapphire library 1",
        "addgene": 233213,
    }
    # Recipes without library identifiers still carry the key, empty.
    assert recipe_metadata("larry")["libraries"] == {}


@pytest.mark.parametrize(
    "libraries",
    [
        (Library("a", "AAAAAA", "A"), Library("a", "CCCCCC", "A again")),
        (Library("a", "AAAAAA", "A"), Library("b", "AAAAAA", "B")),
    ],
)
def test_recipe_rejects_duplicate_libraries(libraries):
    with pytest.raises(ValueError, match="more than once"):
        Recipe(
            name="bad",
            label="Bad",
            description="",
            barcode_length=34,
            stages=(CutadaptStage(adapter="A...B"),),
            libraries=libraries,
        )


def test_recipe_requires_stages_or_an_explanation():
    with pytest.raises(ValueError, match="no cutadapt stages"):
        Recipe(name="bad", label="Bad", description="", barcode_length=8)


def test_register_recipe_refuses_to_clobber():
    with pytest.raises(ValueError, match="already registered"):
        register_recipe(get_recipe("larry"))


def test_no_matching_reads_gives_an_empty_series_not_a_crash():
    """A mismatched recipe must report nothing found, not divide by zero.

    Running the LARRYv2 recipe over LARRY reads is the realistic form of
    this: picking the wrong recipe is the likeliest new failure mode now that
    there is more than one.
    """
    files = {
        "r1": str(HERE / "data/mock_larry_cutadapt_r1.fq"),
        "r2": str(HERE / "data/mock_larry_cutadapt_r2.fq"),
    }
    result = cbu.get_barcodes(files, type="larryv2")

    assert isinstance(result, cbu.CBUclasses.CBUSeries)
    assert len(result) == 0
    assert result.index.names == ["CBC", "Barcode", "UMI"]


def test_thresholds_default_dispr_filter():
    assert Thresholds(min_reads=1, min_umis=1, collapse_distance=1).dispr_filter == 0.1


def test_cli_lists_names():
    """The path a Nextflow config would use to read recipe metadata."""
    out = subprocess.run(
        [sys.executable, "-m", "cbutools", "recipes", "--names-only"],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
        check=True,
    )
    assert out.stdout.split() == RECIPE_NAMES


def test_cli_emits_json():
    out = subprocess.run(
        [sys.executable, "-m", "cbutools", "recipes", "larry"],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
        check=True,
    )
    assert json.loads(out.stdout)["barcode_length"] == 40
