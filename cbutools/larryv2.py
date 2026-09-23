"""LARRYv2 barcode extraction.

LARRYv2 places a 6 nt construct-identifying prefix and a 28 nt barcode motif
between conserved vector sequence. The motif holds 20 random positions separated
by the fixed dinucleotides TG, CA, GT, AG -- fewer random positions and fewer
spacers than LARRY v1's 40-mer, so the two are not interchangeable.

The chemistry is described by the ``larryv2`` entry in :mod:`cbutools.recipes`.
Every construct is read in one pass and the prefix is kept, giving a 34 bp
barcode; demultiplex downstream on its first 6 bases using
``cbutools.recipes.LARRYV2_CONSTRUCTS``.
"""

from .extract import process_recipe
from .recipes import get_recipe


def process_larryv2(
    files,
    valid_CBC=None,
    debug=False,
    cores=1,
    intermediate_dir=None,
):
    """Retrieve CBUSeries from LARRYv2 fastq files

    Processes 10x fastq files (read1 and read2) containing LARRYv2 barcode
    sequences. Cutadapt extracts the sequences matching the barcodes (written
    into temporary files). Then the cell barcodes (CBC), LARRYv2 barcodes and
    UMIs are extracted, counted and assembled into a CBUSeries object.

    Barcodes are 34 bp: the 6 nt construct prefix followed by the 28 nt motif.

    Parameters
    ----------
    files : dict
        Dictionary with paths to fastq files, needs to have 'r1' and 'r2' keys
    valid_CBC : Optional[list, np.ndarray]
        List or array with valid cell barcodes. Supplying this argument is
        recommended as it can dramatically reduce the returned object size.
    debug : bool
        If true the function returns also the temporary directory object, which
        prevents file cleanup and thus allows inspection.
    cores : int
        Number of cores passed to cutadapt.
    intermediate_dir : path-like, optional
        Directory in which cutadapt intermediate files are written. When omitted,
        a unique temporary directory is created in the current working directory.

    Returns
    -------
    CBUSeries or tuple
        Returns a CBUSeries by default. With debug=True, also returns the
        TemporaryDirectory object, or the explicit intermediate directory Path.

    """
    return process_recipe(
        get_recipe("larryv2"),
        files,
        valid_CBC=valid_CBC,
        debug=debug,
        cores=cores,
        intermediate_dir=intermediate_dir,
    )
