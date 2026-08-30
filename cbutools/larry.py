import numpy as np
import pandas as pd
from pathlib import Path
from subprocess import run
import dnaio
from .CBUclasses import CBUSeries
from tempfile import TemporaryDirectory


def process_larry(
    files,
    valid_CBC=None,
    debug=False,
    cores=1,
    intermediate_dir=None,
):
    """Retrieve CBUSeries from LARRY fastq files

    Processes 10x fastq files (read1 and read2) containing LARRY barcode
    sequences. Cutadapt extracts the sequences matching the barcodes (written
    into temporary files). Then the cell barcodes (CBC), LARRY barcodes and UMIs
    are extracted, counted and assembled into a CBUSeries object.

    Parameters
    ----------
    files : dict
        Dictionary with paths to fast files, needs to have 'r1' and 'r2' keys
    valid_CBC : Optional[list, np.ndarray]
        List or array with valid cell barcodes. Supplying this argument is
        recommended as it can dramatically reduce the returned object size.
    debug : bool
        If true the the function returns also the temporary
        directory object, which prevents file cleanup and thus allows inspection.
    cores : int
        Number of cores passed to cutadapt.
    intermediate_dir : path-like, optional
        Directory in which cutadapt intermediate files are written. When omitted,
        a unique temporary directory is created in the current working directory.
        Workflow engines may supply a task-local directory explicitly.

    Returns
    -------
    CBUSeries or tuple
        Returns a CBUSeries by default. With debug=True, also returns the
        TemporaryDirectory object, or the explicit intermediate directory Path.

    """
    temp = None
    if intermediate_dir is None:
        temp = TemporaryDirectory(prefix=".intermediate_cutadapt_", dir="./")
        temp_dir = Path(temp.name).resolve()
    else:
        temp_dir = Path(intermediate_dir).resolve()
        temp_dir.mkdir(parents=True, exist_ok=True)

    cutadapt_threads = ["-j", str(cores)] if cores is not None else []

    # Trimming the files with cutadapt to retain only larry matching barcodes
    command1 = [
        "cutadapt",
        *cutadapt_threads,
        "-g",
        "TTGCTAGGAGAGACCATATG...ATGTCTGGATCCGATATCGC",
        "-o",
        str(temp_dir / "bar1.fq"),
        "-p",
        str(temp_dir / "cbcumi1.fq"),
        str(files["r2"]),
        str(files["r1"]),
        "--discard-untrimmed",
        "--pair-filter=both",
    ]
    # Optionally save the report to json format with  --json=test.cutadapt.json
    out1 = run(command1, capture_output=True, text=True)
    print(out1.stderr)
    print(out1.stdout)
    out1.check_returncode()

    command2 = [
        "cutadapt",
        *cutadapt_threads,
        "-g",
        "NNNNTGNNNNCANNNNACNNNNGANNNNGTNNNNAGNNNN;min_overlap=36",
        "-o",
        str(temp_dir / "bar2.fq"),
        "-p",
        str(temp_dir / "cbcumi2.fq"),
        str(temp_dir / "bar1.fq"),
        str(temp_dir / "cbcumi1.fq"),
        "--discard-untrimmed",
        "--action=retain",
        "--pair-filter=both",
    ]
    out2 = run(command2, capture_output=True, text=True)
    print(out2.stderr)
    print(out2.stdout)
    out2.check_returncode()

    def read_barCBCUMI(x):
        bar = x[0].sequence
        cbcumi = x[1].sequence
        if (len(bar) == 40) and (
            len(cbcumi) == 28
        ):  # Only keeping LARRY barcodes with 40nt and CBC-UMI with 28
            cbc = cbcumi[:16]
            umi = cbcumi[16:28]
            k = (cbc, bar, umi)
            return k
        else:
            return "WRONG_LEN"

    # Accumulating everything in a dictionary with tuples as keys
    barfile = temp_dir / "bar2.fq"
    cbufile = temp_dir / "cbcumi2.fq"
    if valid_CBC is not None:
        counts = {"WRONG_LEN": 0, "CBC_ABSENT": 0}
        with dnaio.open(file1=barfile, file2=cbufile) as reader:
            for record in reader:
                k = read_barCBCUMI(record)
                if k[0] in valid_CBC:
                    if not k in counts:
                        counts[k] = 0
                    counts[k] += 1
                elif k == "WRONG_LEN":
                    counts["WRONG_LEN"] += 1
                else:
                    counts["CBC_ABSENT"] += 1
    else:
        counts = {"WRONG_LEN": 0, "CBC_ABSENT": 0}
        with dnaio.open(file1=barfile, file2=cbufile) as reader:
            for record in reader:
                k = read_barCBCUMI(record)
                if not k in counts:
                    counts[k] = 0
                counts[k] += 1

    wronglen_count = counts["WRONG_LEN"]
    cbcabsent_count = counts["CBC_ABSENT"]
    del counts["WRONG_LEN"]
    del counts["CBC_ABSENT"]

    index = pd.MultiIndex.from_tuples(
        [i for i in counts.keys()], names=("CBC", "Barcode", "UMI")
    )
    counts = CBUSeries([v for v in counts.values()], index=index)

    total_count = counts.sum() + wronglen_count + cbcabsent_count
    valid_count = counts.sum()
    print("============================================================")
    print("Counting LARRY barcodes and CBC UMIs")
    print(
        f"Reads with wrong length: {wronglen_count} ({wronglen_count/total_count*100:.2f}%)"
    )
    print(
        f"Reads with CBC absent in the list: {cbcabsent_count}"
        f"({cbcabsent_count/total_count*100:.2f}%)"
    )
    print(f"Valid reads: {valid_count} ({valid_count/total_count*100:.2f}%)")

    if debug:
        return counts, temp if temp is not None else temp_dir
    else:
        return counts
