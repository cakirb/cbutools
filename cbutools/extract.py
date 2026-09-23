"""Recipe-driven barcode extraction.

One engine, parameterised by a :class:`~cbutools.recipes.Recipe`. It runs the
recipe's cutadapt stages over the read pair, splits each surviving mate into a
cell barcode and a UMI, and counts the resulting ``(CBC, Barcode, UMI)`` triples
into a :class:`~cbutools.CBUclasses.CBUSeries`.

The ``CBUSeries`` is the entire interface to everything downstream. Nothing
here knows what a barcode means.
"""

import pandas as pd
import dnaio
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from .CBUclasses import CBUSeries
from .recipes import CBU_NAMES, Recipe

#: Sentinel counted instead of a triple when a read pair survived cutadapt but
#: has the wrong barcode or CBC-UMI length.
WRONG_LEN = "WRONG_LEN"

#: Sentinel counted instead of a triple when the cell barcode is not in the
#: caller's list of valid cell barcodes.
CBC_ABSENT = "CBC_ABSENT"


def _run_stages(recipe, files, temp_dir, cores):
    """Run the recipe's cutadapt stages and return the final output pair.

    Stage *i* writes ``bar{i}.fq`` and ``cbcumi{i}.fq``. The first stage reads
    the original FASTQs -- read 2 carries the barcode, read 1 the cell barcode
    and UMI -- and each later stage reads its predecessor's output.
    """
    cutadapt_threads = ["-j", str(cores)] if cores is not None else []

    barcode_in = str(files["r2"])
    cbcumi_in = str(files["r1"])

    for number, stage in enumerate(recipe.stages, start=1):
        barcode_out = temp_dir / f"bar{number}.fq"
        cbcumi_out = temp_dir / f"cbcumi{number}.fq"

        command = [
            "cutadapt",
            *cutadapt_threads,
            "-g",
            stage.adapter,
            "-o",
            str(barcode_out),
            "-p",
            str(cbcumi_out),
            barcode_in,
            cbcumi_in,
            "--discard-untrimmed",
            "--pair-filter=both",
        ]
        if stage.action is not None:
            command.append(f"--action={stage.action}")

        # Optionally save the report to json format with --json=test.cutadapt.json
        out = run(command, capture_output=True, text=True)
        print(out.stderr)
        print(out.stdout)
        out.check_returncode()

        barcode_in = str(barcode_out)
        cbcumi_in = str(cbcumi_out)

    return Path(barcode_in), Path(cbcumi_in)


def _count_triples(recipe, barfile, cbufile, valid_CBC):
    """Count (CBC, Barcode, UMI) triples, plus the two rejection counters."""

    def read_barCBCUMI(x):
        bar = x[0].sequence
        cbcumi = x[1].sequence
        if recipe.accepts_barcode(bar) and (len(cbcumi) == recipe.cbcumi_length):
            cbc, umi = recipe.split_cbcumi(cbcumi)
            return (cbc, bar, umi)
        else:
            return WRONG_LEN

    counts = {WRONG_LEN: 0, CBC_ABSENT: 0}
    with dnaio.open(file1=barfile, file2=cbufile) as reader:
        if valid_CBC is not None:
            for record in reader:
                k = read_barCBCUMI(record)
                if k == WRONG_LEN:
                    counts[WRONG_LEN] += 1
                elif k[0] in valid_CBC:
                    if not k in counts:
                        counts[k] = 0
                    counts[k] += 1
                else:
                    counts[CBC_ABSENT] += 1
        else:
            for record in reader:
                k = read_barCBCUMI(record)
                if not k in counts:
                    counts[k] = 0
                counts[k] += 1

    wronglen_count = counts.pop(WRONG_LEN)
    cbcabsent_count = counts.pop(CBC_ABSENT)
    return counts, wronglen_count, cbcabsent_count


def _to_series(counts):
    """Build a CBUSeries from the counted triples, empty input included."""
    if not counts:
        # from_tuples cannot infer a level count from an empty list, and an
        # empty result is exactly what a mismatched recipe produces -- the case
        # where a clear answer matters most.
        index = pd.MultiIndex.from_arrays([[], [], []], names=CBU_NAMES)
        return CBUSeries([], index=index, dtype="int64")

    index = pd.MultiIndex.from_tuples([i for i in counts.keys()], names=CBU_NAMES)
    return CBUSeries([v for v in counts.values()], index=index)


def _report(recipe, counts, wronglen_count, cbcabsent_count):
    """Print the rejection breakdown. The first thing anyone looks at."""
    valid_count = int(counts.sum())
    total_count = valid_count + wronglen_count + cbcabsent_count

    def pct(n):
        return f"{n/total_count*100:.2f}%" if total_count else "n/a"

    print("============================================================")
    print(f"Counting {recipe.label} barcodes and CBC UMIs")
    print(f"Reads with wrong length: {wronglen_count} ({pct(wronglen_count)})")
    print(
        f"Reads with CBC absent in the list: {cbcabsent_count}"
        f"({pct(cbcabsent_count)})"
    )
    print(f"Valid reads: {valid_count} ({pct(valid_count)})")

    if valid_count == 0:
        print(
            f"WARNING: no valid reads for recipe '{recipe.name}'. Check that the "
            f"reads really are {recipe.label} ({recipe.barcode_length} bp barcode, "
            f"{recipe.chemistry}) before trusting an empty result."
        )


def process_recipe(
    recipe,
    files,
    valid_CBC=None,
    debug=False,
    cores=1,
    intermediate_dir=None,
):
    """Retrieve a CBUSeries from a fastq pair using a barcode recipe

    Processes 10x fastq files (read1 and read2) containing lineage barcodes.
    Cutadapt extracts the sequences matching the recipe's adapters (written into
    temporary files). Then the cell barcodes (CBC), lineage barcodes and UMIs are
    extracted, counted and assembled into a CBUSeries object.

    Parameters
    ----------
    recipe : Recipe or str
        The recipe to apply, or the name of a registered one.
    files : dict
        Dictionary with paths to fastq files, needs to have 'r1' and 'r2' keys
    valid_CBC : Optional[list, np.ndarray]
        List or array with valid cell barcodes. Supplying this argument is
        recommended as it can dramatically reduce the returned object size.
    debug : bool
        If true the function returns also the temporary directory object, which
        prevents file cleanup and thus allows inspection.
    cores : int
        Number of cores passed to cutadapt. Pass the number your scheduler
        actually reserved; the default of 1 does not mean "all cores".
    intermediate_dir : path-like, optional
        Directory in which cutadapt intermediate files are written. When omitted,
        a unique temporary directory is created in the current working directory.
        Workflow engines may supply a task-local directory explicitly.

    Raises
    ------
    NotImplementedError
        If the recipe is registered for its metadata but not yet runnable.

    Returns
    -------
    CBUSeries or tuple
        Returns a CBUSeries by default. With debug=True, also returns the
        TemporaryDirectory object, or the explicit intermediate directory Path.

    """
    if not isinstance(recipe, Recipe):
        from .recipes import get_recipe

        recipe = get_recipe(recipe)

    if not recipe.runnable:
        raise NotImplementedError(
            f"Recipe {recipe.name!r} is not runnable yet. {recipe.pending}"
        )

    temp = None
    if intermediate_dir is None:
        temp = TemporaryDirectory(prefix=".intermediate_cutadapt_", dir="./")
        temp_dir = Path(temp.name).resolve()
    else:
        temp_dir = Path(intermediate_dir).resolve()
        temp_dir.mkdir(parents=True, exist_ok=True)

    barfile, cbufile = _run_stages(recipe, files, temp_dir, cores)
    counts, wronglen_count, cbcabsent_count = _count_triples(
        recipe, barfile, cbufile, valid_CBC
    )
    counts = _to_series(counts)
    _report(recipe, counts, wronglen_count, cbcabsent_count)

    if debug:
        return counts, temp if temp is not None else temp_dir
    else:
        return counts
