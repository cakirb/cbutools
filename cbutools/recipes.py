"""Barcode recipes.

A recipe is the complete description of one lineage-barcoding chemistry. It has
two halves, and both matter:

* **Extraction** -- the cutadapt stages that cut the barcode out of a read pair,
  and the lengths used to split the mate into a cell barcode and a UMI. This
  half is consumed by :mod:`cbutools.extract`.
* **Metadata** -- barcode length, sensible starting thresholds for filtering
  and Hamming collapse, and whether a
  whitelist of expected barcodes exists. This half is consumed by whatever drives
  ``cbutools``; a workflow engine cannot infer any of it from the returned
  ``CBUSeries``, and getting it wrong is silent rather than loud. A recipe is not
  fully described by its extraction logic.

Recipes are data. Adding a chemistry should mean adding a :class:`Recipe` to
``_RECIPES`` below, not writing new code.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Tuple

#: Near-identical barcodes differ by substitution only. Valid because every
#: recipe's barcodes have one fixed length. This is what
#: ``CBUSeries.filter_by_hamming`` implements, and it is the only collapse metric
#: cbutools supports; recipe metadata reports it as ``collapse_metric``.
HAMMING = "hamming"

#: Index level names of a ``CBUSeries``. The one fixed part of the contract.
CBU_NAMES = ("CBC", "Barcode", "UMI")


@dataclass(frozen=True)
class CutadaptStage:
    """One cutadapt pass over the (barcode, CBC-UMI) read pair.

    Stages run in order. The first reads the original FASTQ pair; each later one
    reads the previous stage's output. Every stage runs with
    ``--discard-untrimmed --pair-filter=both``, so a read pair survives only
    while it keeps matching.

    Parameters
    ----------
    adapter : str
        Value passed to cutadapt's ``-g``, verbatim. Adapter parameters such as
        ``;min_overlap=36`` and ``;required`` belong here.
    action : str, optional
        Value passed to cutadapt's ``--action``. ``None`` leaves the default
        (trim), which removes the matched region; ``"retain"`` keeps it.
    description : str
        What this stage is for. Shown in ``--help``-style metadata dumps.
    """

    adapter: str
    action: Optional[str] = None
    description: str = ""


@dataclass(frozen=True)
class Thresholds:
    """Filtering thresholds for one recipe under one stringency profile.

    ``collapse_distance`` is a Hamming distance over the whole barcode.

    **It is an exclusive threshold**, matching ``hamming_filter(min_distance=)``
    and nf-larry's ``min_hamming``: sequences strictly *closer* than this are
    collapsed. So ``collapse_distance=3`` acts on pairs at d<=2, and a published
    "collapse within Hamming distance 4" becomes ``collapse_distance=5``. This
    off-by-one is the easiest thing in the package to get wrong, and getting it
    wrong is silent.
    """

    min_reads: int
    min_umis: int
    collapse_distance: int
    dispr_filter: float = 0.1


@dataclass(frozen=True)
class Library:
    """One officially distributed barcode library of a recipe.

    Some chemistries ship as several libraries whose barcodes start with a
    short library identifier. The identifier stays inside the extracted barcode;
    this record only names it, so a caller can say which library a barcode came
    from.

    Parameters
    ----------
    key : str
        Short machine-friendly name, e.g. ``tsapphire1``.
    identifier : str
        The identifier sequence at the start of the barcode.
    label : str
        The distributor's name for the library.
    addgene : int, optional
        Addgene catalogue number of the pooled library.
    """

    key: str
    identifier: str
    label: str
    addgene: Optional[int] = None


@dataclass(frozen=True)
class Recipe:
    """A complete barcoding chemistry description.

    Parameters
    ----------
    name : str
        Recipe selector, as passed to ``get_barcodes(type=...)``.
    label : str
        Human-readable name, used in printed extraction summaries.
    description : str
        One-line summary of the chemistry.
    barcode_length : int
        Expected barcode length in bases. Every barcode has this length, which
        is what makes Hamming collapse valid.
    stages : tuple of CutadaptStage
        The extraction passes. Empty only for a recipe that is not yet runnable,
        in which case ``pending`` says what is missing.
    cbc_length, umi_length : int
        How the mate read is split. The defaults are 10x Genomics 3' v3/v3.1
        (16 bp cell barcode, 12 bp UMI) and are a property of the sequencing
        chemistry, not of the barcoding tool, so they carry across recipes.
    chemistry : str
        The single-cell chemistry the recipe is valid for.
    defaults, lenient : Thresholds
        Starting thresholds for a stringent and a permissive pass. Thresholds do
        not transfer between recipes -- a collapse distance of 3 is mild on a
        40 bp barcode and would be drastic on a short one.
    libraries : tuple of Library
        Officially distributed libraries, for recipes whose barcodes start with
        a library identifier. Empty otherwise. Barcodes with an identifier not
        listed here still extract; the list names libraries, it does not filter.
    whitelist : str, optional
        Where the set of expected barcodes comes from, if one exists. A
        description or URL, not a bundled file -- ``cbutools`` ships no
        whitelists.
    reference : str
        Citation for the chemistry.
    notes : tuple of str
        Caveats a caller should read before trusting the recipe.
    pending : str, optional
        If set, the recipe is registered for its metadata but cannot be run, and
        this says exactly what is missing. Attempting extraction raises
        ``NotImplementedError`` carrying this text.
    """

    name: str
    label: str
    description: str
    barcode_length: int
    stages: Tuple[CutadaptStage, ...] = ()
    cbc_length: int = 16
    umi_length: int = 12
    chemistry: str = "10x Genomics 3' v3 / v3.1"
    defaults: Thresholds = field(
        default=Thresholds(min_reads=2, min_umis=2, collapse_distance=2)
    )
    lenient: Thresholds = field(
        default=Thresholds(min_reads=1, min_umis=1, collapse_distance=2)
    )
    libraries: Tuple[Library, ...] = ()
    whitelist: Optional[str] = None
    reference: str = ""
    notes: Tuple[str, ...] = ()
    pending: Optional[str] = None

    def __post_init__(self):
        if not isinstance(self.barcode_length, int) or self.barcode_length <= 0:
            raise ValueError(
                f"Recipe {self.name!r} needs a positive fixed barcode_length; "
                "Hamming collapse is undefined for sequences of unequal length"
            )
        for attr in ("key", "identifier"):
            values = [getattr(library, attr) for library in self.libraries]
            if len(set(values)) != len(values):
                raise ValueError(
                    f"Recipe {self.name!r} lists a library {attr} more than once"
                )
        if not self.stages and not self.pending:
            raise ValueError(
                f"Recipe {self.name!r} has no cutadapt stages and no 'pending' "
                "explanation; one or the other is required"
            )

    @property
    def cbcumi_length(self) -> int:
        """Expected length of the mate read carrying the cell barcode and UMI."""
        return self.cbc_length + self.umi_length

    @property
    def runnable(self) -> bool:
        """Whether extraction can actually be performed with this recipe."""
        return self.pending is None

    def split_cbcumi(self, cbcumi: str) -> Tuple[str, str]:
        """Split a mate read into (cell barcode, UMI)."""
        return cbcumi[: self.cbc_length], cbcumi[self.cbc_length : self.cbcumi_length]

    def accepts_barcode(self, barcode: str) -> bool:
        """Whether an extracted barcode has an acceptable length."""
        return len(barcode) == self.barcode_length

    def to_dict(self) -> Dict:
        """Return a JSON-serialisable description of the recipe.

        This is the form a workflow engine should consume: it carries the
        metadata needed to configure downstream collapse and filtering without
        importing ``cbutools``.
        """
        return {
            "name": self.name,
            "label": self.label,
            "description": self.description,
            "chemistry": self.chemistry,
            "barcode_length": self.barcode_length,
            "cbc_length": self.cbc_length,
            "umi_length": self.umi_length,
            "collapse_metric": HAMMING,
            "defaults": asdict(self.defaults),
            "lenient": asdict(self.lenient),
            "libraries": {
                library.key: {
                    "identifier": library.identifier,
                    "label": library.label,
                    "addgene": library.addgene,
                }
                for library in self.libraries
            },
            "whitelist": self.whitelist,
            "reference": self.reference,
            "notes": list(self.notes),
            "runnable": self.runnable,
            "pending": self.pending,
            "stages": [
                {
                    "adapter": stage.adapter,
                    "action": stage.action,
                    "description": stage.description,
                }
                for stage in self.stages
            ],
        }


# --------------------------------------------------------------------------- #
# The recipes themselves.
# --------------------------------------------------------------------------- #

_LARRY = Recipe(
    name="larry",
    label="LARRY",
    description=(
        "LARRY lineage barcode, 40 bp, expressed downstream of the GFP coding "
        "sequence and immediately before the polyadenylation signal."
    ),
    barcode_length=40,
    stages=(
        CutadaptStage(
            adapter="TTGCTAGGAGAGACCATATG...ATGTCTGGATCCGATATCGC",
            description="Cut the barcode out from between its flanking constant regions.",
        ),
        CutadaptStage(
            adapter="NNNNTGNNNNCANNNNACNNNNGANNNNGTNNNNAGNNNN;min_overlap=36",
            action="retain",
            description=(
                "Require the barcode motif itself: random 4-mers separated by "
                "fixed dinucleotide spacers. Retained rather than trimmed."
            ),
        ),
    ),
    # Values in use in the nf-larry pipeline, tuned on real data. They are not
    # the same as the defaults baked into get_barcodes_report's signature.
    defaults=Thresholds(min_reads=8, min_umis=3, collapse_distance=3, dispr_filter=0.1),
    lenient=Thresholds(min_reads=1, min_umis=2, collapse_distance=3, dispr_filter=0.1),
    whitelist="Library-specific; nf-larry ships one as data/larry_whitelist.csv.",
    reference="Weinreb et al., Science 2020, doi:10.1126/science.aaw3381",
)


# LARRYv2 read structure on 10x 3' v3, read 2 (cloneRanger's own config states
# read_feature_bc: R2, read_cellular_bc: R1, matching LARRY v1):
#
#   [2-4 nt stagger] GCTAGGAGAGACCATATGGGATC CGAT TCCAGT NNNNTGNNNNCANNNNGTNNNNAGNNNN ATCGCCACC...
#                    |<---- anchor, 23 ---->|      |<-6->| |<------ motif, 28 ------->| |<- flank ->|
#
# The stagger is 2-4 nt, measured over 300k reads: the anchor starts at offset
# 2/3/4 in a near-even three-way split, and TCCAGT at anchor_start + 27 in every
# bin. Both extract_larryv2.py's docstring and REQ-72680's stage-1 report say
# "0-2 nt"; they are off by two. It makes no difference to a non-anchored
# cutadapt match, but it would if the 5' end were ever anchored.
#
# The anchor is the tail of the paper's Trueseq_LARRY nested-PCR forward primer,
# and overlaps LARRY v1's own 5' flank (GCTAGGAGAGACCATATG) -- this is conserved
# vector sequence. The 6 nt after CGAT is the library identifier: it says which
# barcode library the barcode came from; see LARRYV2_LIBRARIES.
LARRYV2_ANCHOR = "GCTAGGAGAGACCATATGGGATCCGAT"

# First 9 nt of the constant region that follows the barcode. Measured over
# 287k pattern-matching reads, the exact-match rate is 99.02% at 6 nt, 98.91% at
# 9 and 98.53% at 15, so the region is genuinely invariant and the extra length
# is nearly free; the ~1% residual is ordinary substitution and indel error, not
# a second library. Nine keeps headroom on shorter read 2 lengths -- extend to
# ATCGCCACCGTG if you want more specificity and your reads allow it.
LARRYV2_FLANK3 = "ATCGCCACC"

# 20 random positions with four fixed dinucleotides. Note the spacers differ
# from LARRY v1's: TG, CA, GT, AG here, against TG, CA, AC, GA, GT, AG there.
LARRYV2_MOTIF = "NNNNTGNNNNCANNNNGTNNNNAGNNNN"

#: The official LARRYv2 barcode libraries distributed by Addgene
#: (addgene.org/pooled-library/rodriguez-fraticelli-larryv2-barcoding/). Each
#: barcode starts with its library's 6 nt identifier; the depositor's figure
#: gives them as "GTTCCA-barcode for LARRYv2 T-Sapphire library 1" and so on.
#: The two libraries of one fluorophore share a backbone and differ only in the
#: barcode oligo pool. mScarlet (233697) was listed as back-ordered when checked
#: on 2026-11-08.
#:
#: These are not separate recipes. Every library is read in one pass and the
#: identifier stays inside the barcode, so the library is part of the clone
#: identity and a multiplexed sample comes back whole. Read the library from
#: ``barcode[:6]``.
#:
#: Identifiers not listed here extract just the same. Singh et al. 2025 (STRACK),
#: for example, also used the lab's unpublished GCTGAA and CTAGAT libraries;
#: callers should report such identifiers by sequence rather than guess a name.
LARRYV2_LIBRARIES: Dict[str, Library] = {
    library.key: library
    for library in (
        Library("tsapphire1", "GTTCCA", "LARRYv2-T-Sapphire library 1", 233213),
        Library("tsapphire2", "TGCTAA", "LARRYv2-T-Sapphire library 2", 233214),
        Library("egfp1", "TCCAGT", "LARRYv2-EGFP library 1", 233215),
        Library("egfp2", "TGATTG", "LARRYv2-EGFP library 2", 233216),
        Library("mscarlet", "AAGATT", "LARRYv2-mScarlet library", 233697),
    )
}

#: Minimum pairwise Hamming distance between the library identifiers: 3, between
#: TGCTAA (tsapphire2) and TGATTG (egfp2). Collapse runs on the whole 34 bp
#: string, identifier included, so it cannot fuse barcodes from two libraries
#: while collapse_distance <= this value. tests/test_recipes.py asserts that, so
#: adding a library with a closer identifier fails there rather than in clone
#: calls.
LARRYV2_MIN_PREFIX_DISTANCE = 3

_LARRYV2 = Recipe(
    name="larryv2",
    label="LARRYv2",
    description=(
        "LARRYv2 lineage barcode, 34 bp: the 6 nt library identifier followed "
        "by the 28 nt barcode motif, expressed from an "
        "EF1a-Reporter-WPRE-barcode-pA construct. Every library is read in "
        "one pass and the identifier is kept, so a multiplexed sample comes "
        "back whole."
    ),
    barcode_length=34,
    stages=(
        CutadaptStage(
            adapter=(
                f"{LARRYV2_ANCHOR};min_overlap={len(LARRYV2_ANCHOR)};required"
                f"...{LARRYV2_FLANK3};min_overlap={len(LARRYV2_FLANK3)};required"
            ),
            description=(
                "Cut the library identifier and barcode out from between the "
                "vector anchor and the downstream constant region."
            ),
        ),
        CutadaptStage(
            # e=0 so the eight fixed spacer bases must match exactly. At
            # cutadapt's default 10% rate a 34 nt adapter would tolerate three
            # errors, which is more than half the fixed content of the motif.
            # The identifier is left as NNNNNN, so a library outside
            # LARRYV2_LIBRARIES still extracts.
            adapter=f"{'N' * 6}{LARRYV2_MOTIF};min_overlap=34;e=0",
            action="retain",
            description=(
                "Require the barcode motif: six unconstrained identifier bases "
                "then 20 random positions separated by the fixed dinucleotides "
                "TG, CA, GT, AG. Retained rather than trimmed."
            ),
        ),
    ),
    # min_reads and min_umis are cloneRanger's reads_cutoff: 5 and umi_cutoff: 3.
    # collapse_distance is NOT cloneRanger's: its hamming_distance: 4 (d<=4) is
    # safe there because it keeps the library as a label beside the barcode.
    # Here the identifier is inside the compared string, and at d<=4 two
    # barcodes from libraries 3 apart whose random parts differ by one base
    # would fuse.
    # 3 (d<=2) makes that impossible, and still absorbs essentially all real
    # error: at REQ-72680's ~2.6e-3 per base, about 5% of barcodes carry one
    # error, 0.13% two, and 0.002% three or more.
    defaults=Thresholds(min_reads=5, min_umis=3, collapse_distance=3, dispr_filter=0.1),
    lenient=Thresholds(min_reads=1, min_umis=2, collapse_distance=3, dispr_filter=0.1),
    libraries=tuple(LARRYV2_LIBRARIES.values()),
    whitelist=(
        "Library-specific; supplied with the Addgene pooled library. Whitelist "
        "entries should be the 34 bp identifier+motif, so one list spans libraries."
    ),
    reference=(
        "Rodriguez-Fraticelli lab LARRYv2 libraries; library names and "
        "identifiers from Addgene "
        "(addgene.org/pooled-library/rodriguez-fraticelli-larryv2-barcoding/). "
        "Read structure cross-checked against cloneRanger "
        "config/larry_config.yaml (github.com/dfernandezperez/cloneRanger)."
    ),
    notes=(
        "The first 6 bases of each barcode identify the Addgene library "
        "(see libraries). The same random bases under two different "
        "identifiers are two clones, which is why the identifier is kept.",
        "Identifiers outside the listed libraries extract too, e.g. the "
        "unpublished GCTGAA and CTAGAT libraries in Singh et al. 2025. Report "
        "them by sequence; the list names libraries, it does not filter them.",
        "collapse_distance is 3, not cloneRanger's d<=4 (which would be 5 "
        "here). Collapse runs on the whole string, identifier included, and the "
        "closest identifiers (TGCTAA tsapphire2, TGATTG egfp2) are only 3 apart. "
        "If the caller ever groups barcodes by identifier before collapsing -- "
        "as cloneRanger does -- 5 becomes safe within each group.",
        "The 20 random positions are not equimolar (roughly A 24%, C 19%, G 23%, "
        "T 34%), so effective diversity is about 6.7e11 rather than 4^20 = 1.1e12.",
        "These thresholds are for single-cell scale (~1e3-1e5 barcodes). Do not "
        "reuse them for plasmid library diversity estimation, where the "
        "library is three orders of magnitude larger.",
        "The flanks were verified on plasmid amplicon data (REQ-72680, EGFP "
        "library 1) and against cloneRanger, not on a 10x read pair. The read "
        "layout follows cloneRanger's read_feature_bc: R2 / read_cellular_bc: "
        "R1; confirm against a real 10x LARRYv2 sample before a production run.",
    ),
)


_RECIPES: Dict[str, Recipe] = {
    recipe.name: recipe
    for recipe in (
        _LARRY,
        _LARRYV2,
    )
}


def list_recipes(include_pending: bool = False) -> List[str]:
    """List registered recipe names.

    Parameters
    ----------
    include_pending : bool
        Whether to include recipes registered for their metadata but not yet
        runnable. Off by default, so that a caller offering the list as a choice
        does not offer something that will fail.

    Returns
    -------
    list of str
    """
    return sorted(
        name for name, recipe in _RECIPES.items() if include_pending or recipe.runnable
    )


def get_recipe(name: str) -> Recipe:
    """Look up a recipe by name.

    Raises
    ------
    KeyError
        If no such recipe is registered. The message lists what is available.
    """
    try:
        return _RECIPES[name]
    except KeyError:
        raise KeyError(
            f"Unknown barcode recipe {name!r}. Available: "
            f"{', '.join(list_recipes(include_pending=True))}"
        ) from None


def recipe_metadata(name: Optional[str] = None) -> Dict:
    """Return JSON-serialisable recipe metadata.

    Parameters
    ----------
    name : str, optional
        A single recipe. When omitted, every registered recipe is returned,
        keyed by name.

    Returns
    -------
    dict
    """
    if name is not None:
        return get_recipe(name).to_dict()
    return {
        recipe_name: recipe.to_dict() for recipe_name, recipe in sorted(_RECIPES.items())
    }


def register_recipe(recipe: Recipe, overwrite: bool = False) -> Recipe:
    """Register a recipe at runtime.

    Mostly useful for trying a chemistry out before adding it to this module.

    Raises
    ------
    ValueError
        If the name is taken and ``overwrite`` is false.
    """
    if recipe.name in _RECIPES and not overwrite:
        raise ValueError(
            f"Recipe {recipe.name!r} is already registered; pass overwrite=True "
            "to replace it"
        )
    _RECIPES[recipe.name] = recipe
    return recipe


def _main(argv=None):
    """Print recipe metadata as JSON, for consumers that are not Python.

    Reached through ``python -m cbutools recipes``; see ``cbutools/__main__.py``.
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m cbutools recipes",
        description="Print cbutools barcode recipe metadata as JSON.",
    )
    parser.add_argument(
        "name",
        nargs="?",
        help="Recipe to describe. Omit to describe all of them.",
    )
    parser.add_argument(
        "--names-only",
        action="store_true",
        help="Print just the runnable recipe names, one per line.",
    )
    parser.add_argument(
        "--include-pending",
        action="store_true",
        help="With --names-only, also list recipes that are not yet runnable.",
    )
    args = parser.parse_args(argv)

    if args.names_only:
        for name in list_recipes(include_pending=args.include_pending):
            print(name)
        return 0

    print(json.dumps(recipe_metadata(args.name), indent=2))
    return 0

