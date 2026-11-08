# cbutools

cbutools recovers expressed lineage barcodes from 10x Genomics scRNA-seq data and assigns them to single cells.

It was developed by Iwo Kucinski in collaboration with the Berthold Göttgens lab at the Stem Cell Institute, University of Cambridge.

## Installation

```bash
git clone https://github.com/Iwo-K/cbutools
cd cbutools
pip install .
```

Optional extras:

```bash
pip install '.[report]'   # jupyter, papermill and nbconvert, for .ipynb/.html reports
pip install '.[dev]'      # the report extras plus pytest and black
```

`cutadapt` is a required dependency and must be on `PATH`.

## How to use

Call `get_barcodes` with a dictionary of FASTQ paths (keys `'r1'` and `'r2'`) and the recipe for your barcoding chemistry:

```python
import cbutools as cbu

counts = cbu.get_barcodes({"r1": read1, "r2": read2}, type="larry")
```

The result is a `CBUSeries` holding read counts for each CBC-Barcode-UMI combination. It can then be filtered and used to assign barcodes to cells. For a worked example, see the report template at `cbutools/templates/larry1_template.ipynb`.

You can also pass a list of valid cell barcodes, for example the real cells (not empty droplets) from a Cell Ranger run. This is recommended, because it can shrink the returned object a lot:

```python
counts = cbu.get_barcodes(files, type="larry", valid_CBC=cells)
```

Workflow engines can pass their allocated cores and a task-local directory for intermediate files:

```python
cbu.get_barcodes(
    files={"r1": read1, "r2": read2},
    type="larry",
    cores=task_cpus,
    intermediate_dir=".intermediate_cutadapt",
)
```

Pass the number of cores your scheduler actually reserved. The default of 1 does not mean "all cores". If `intermediate_dir` is omitted, a unique temporary directory is used and cleaned up automatically.

### CBUSeries

This behaves like a pandas `Series` with a `MultiIndex` of `(CBC, Barcode, UMI)`, and adds these methods:

- `plot_barcode_no` plots the number of barcodes observed per CBC.
- `plot_hist` plots a histogram of read counts per CBC, Barcode, UMI, or a combination of them.
- `filter_by_reads` keeps entries with at least a minimum number of reads per CBC, Barcode, UMI, or a combination of them.
- `filter_by_hamming` collapses Barcodes (or CBCs, or UMIs) that are closer than a minimum Hamming distance, keeping the more abundant one.
- `summary` calculates and prints basic barcode statistics.
- `save_barcodes` saves to a CSV file that `load_barcodes` can read back.
- `count_UMI` counts UMIs per CBC-Barcode combination and returns a `CBSeries`.

### CBSeries

This behaves like a pandas `Series` with a `MultiIndex` of `(CBC, Barcode)`, and adds these methods:

- `plot_barcode_no` plots the number of barcodes observed per CBC.
- `plot_hist` plots a histogram of UMI counts per CBC, Barcode, or both.
- `filter_by_UMI` keeps entries with at least a minimum number of UMIs per CBC, Barcode, or both.
- `save_barcodes` saves to a CSV file that `load_barcodes` can read back.
- `summary` calculates and prints basic barcode statistics.
- `assign_barcodes` returns a pandas `DataFrame` of the barcodes assigned to each cell.

## Recipes

A recipe describes one barcoding chemistry. Choose one with the `type` argument:

```python
cbu.list_recipes()   # ['larry', 'larryv2']
cbu.get_barcodes(files, type="larryv2")
```

| recipe | barcode | reference |
|---|---|---|
| `larry` | 40 bp | Weinreb et al., *Science* 2020 |
| `larryv2` | 34 bp | Rodriguez-Fraticelli lab LARRYv2 libraries |

Every recipe has fixed-length barcodes, so near-identical barcodes are collapsed by Hamming distance (`CBUSeries.filter_by_hamming`).

All recipes currently assume 10x Genomics 3′ v3 / v3.1 chemistry: read 2 carries the barcode, and read 1 carries the 16 bp cell barcode followed by the 12 bp UMI.

### A recipe is more than its extraction logic

Knowing how to cut a barcode out of a read isn't enough to process it correctly afterwards. You also need to know:

- how long the barcode is,
- which thresholds are sensible.

None of that can be inferred from the returned `CBUSeries`, so each recipe carries it alongside its cutadapt stages:

```python
r = cbu.get_recipe("larryv2")
r.barcode_length     # 34
r.defaults           # Thresholds(min_reads=5, min_umis=3, collapse_distance=3, dispr_filter=0.1)
r.lenient            # the permissive profile
```

**Thresholds don't transfer between recipes.** A collapse distance means different things on barcodes of different lengths and designs, so each recipe declares its own stringent `defaults` and permissive `lenient` profile.

**`collapse_distance` is exclusive.** Barcodes strictly closer than it are collapsed, so `collapse_distance=3` acts on pairs at d≤2. This is the same convention as `filter_by_hamming(min_distance=...)`. Under it, a published "collapse within Hamming distance 4" becomes `collapse_distance=5`. This is the easiest thing in the package to get wrong, and a wrong value fails silently.

Tools that aren't Python can read the same metadata as JSON:

```bash
python -m cbutools recipes                 # every recipe
python -m cbutools recipes larry           # one recipe
python -m cbutools recipes --names-only    # runnable recipe names, one per line
```

### LARRY (`type="larry"`)

LARRY barcodes were introduced by [Weinreb et al., *Science* 2020](https://doi.org/10.1126/science.aaw3381) for lineage tracing of haematopoietic stem and progenitor cells. The barcode sits in the mRNA downstream of the GFP coding sequence, just before the polyadenylation signal.

It is flanked by `TTGCTAGGAGAGACCATATG...ATGTCTGGATCCGATATCGC` and has the form:

```
NNNNTGNNNNCANNNNACNNNNGANNNNGTNNNNAGNNNN
```

The barcode is a fixed 40 bp, so Hamming collapse applies. `process_larry` remains available as a named entry point.

### LARRYv2 (`type="larryv2"`)

LARRYv2 is the [Rodriguez-Fraticelli lab](https://www.addgene.org/pooled-library/rodriguez-fraticelli-larryv2-barcoding/) successor to LARRY. Read 2 has this structure:

```
[2-4 nt stagger] GCTAGGAGAGACCATATGGGATC CGAT TCCAGT NNNNTGNNNNCANNNNGTNNNNAGNNNN ATCGCCACC...
                 |<---- anchor, 23 ---->|      |<-6->| |<------ motif, 28 ------->| |<- flank ->|
```

The motif has 20 random positions and four fixed dinucleotides. LARRY v1's 40-mer has 28 random positions and six, so the two are not interchangeable.

The 6 nt after `CGAT` is the library identifier: it says which barcode library the barcode came from. The official libraries are those [distributed by Addgene](https://www.addgene.org/pooled-library/rodriguez-fraticelli-larryv2-barcoding/) (`cbutools.recipes.LARRYV2_LIBRARIES`):

| key | Addgene library | Addgene ID | identifier |
|---|---|---|---|
| `tsapphire1` | LARRYv2-T-Sapphire library 1 | 233213 | `GTTCCA` |
| `tsapphire2` | LARRYv2-T-Sapphire library 2 | 233214 | `TGCTAA` |
| `egfp1` | LARRYv2-EGFP library 1 | 233215 | `TCCAGT` |
| `egfp2` | LARRYv2-EGFP library 2 | 233216 | `TGATTG` |
| `mscarlet` | LARRYv2-mScarlet library | 233697 | `AAGATT` |

The same table is in the recipe metadata, under `libraries` (`python -m cbutools recipes larryv2`).

Barcodes come back as 34 bp: the library identifier followed by the motif. Every library is read in one pass, as the authors' own pipeline (cloneRanger) does. The identifier is kept because it is part of the clone identity: the same random bases under two libraries are two different clones. To split a multiplexed sample by library, use the first 6 bases of the barcode.

Barcodes with an identifier that isn't in the table extract too, for example the unpublished `GCTGAA` and `CTAGAT` libraries used by Singh et al. 2025. Report them by their sequence.

**Thresholds.** `min_reads` and `min_umis` default to cloneRanger's `reads_cutoff: 5` and `umi_cutoff: 3`.

`collapse_distance` is 3 (merge at d≤2) rather than cloneRanger's d≤4, and this is deliberate:

- cloneRanger keeps the library as a separate label, so its collapse never compares identifiers.
- Here the identifier is inside the compared string, and the closest pair of identifiers (`TGCTAA` tsapphire2, `TGATTG` egfp2) is only 3 apart.
- So at d≤4, two barcodes from different libraries whose random bases differ by one would be fused into one clone. At d≤2 that can't happen.

The cost is negligible. At the ~2.6×10⁻³ per-base error rate measured on real LARRYv2 reads, about 5% of barcodes carry one error and 0.13% carry two, and both are absorbed at d≤2. Barcodes with three or more errors are about 2 in 100,000.

Caveats:

- The flanks were verified on plasmid amplicon data and against cloneRanger's config, not on a 10x read pair. The read layout follows cloneRanger's `read_feature_bc: R2` / `read_cellular_bc: R1`. Confirm it against a real 10x LARRYv2 sample before a production run.

### Adding a recipe

Recipes are data. To add a chemistry, add a `Recipe` to `_RECIPES` in `cbutools/recipes.py`. You don't need to write new code, because the extraction engine in `cbutools/extract.py` is shared.

A `Recipe` lists:

- its cutadapt stages. Each runs with `--discard-untrimmed --pair-filter=both`. The first reads the original FASTQ pair, and each later one reads the previous stage's output.
- the barcode, CBC and UMI lengths.
- the `defaults` and `lenient` threshold profiles.

`Recipe` rejects a recipe without a fixed barcode length, because Hamming collapse needs one. To try a chemistry out before adding it to the module, register it at runtime with `register_recipe`.

## Reports

cbutools can generate reports with jupyter and papermill (install the `[report]` or `[dev]` extras). Templates are stored in `cbutools/templates/` and run with `get_barcodes_report`. To add a template, put an `.ipynb` notebook in that directory and name it in `get_barcodes_report`.

## Testing

```bash
pytest
```

The extraction tests need `cutadapt` on `PATH`. When you're testing reports, setting pytest's temporary directory makes the generated reports easier to inspect:

```bash
pytest --basetemp=mydir
```
