"""LARRY barcode extraction.

The chemistry itself is described by the ``larry`` entry in
:mod:`cbutools.recipes`; the extraction is done by
:func:`cbutools.extract.process_recipe`. This module remains as the named entry
point that existing callers import.
"""

from .extract import process_recipe
from .recipes import get_recipe


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
    return process_recipe(
        get_recipe("larry"),
        files,
        valid_CBC=valid_CBC,
        debug=debug,
        cores=cores,
        intermediate_dir=intermediate_dir,
    )
