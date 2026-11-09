# Python_local

A small library for making PSI-BLAST position-specific scoring matrices (PSSMs) from protein FASTA files,
to use as features for protein classifiers (for example peroxisome targeting; see the `Machine-Learning`
repository).

**Status:** work in progress. `pssm.py` is the current version; `deneme.py` and `pssm_class.ipynb` are
earlier drafts (*deneme* = "test").

## Files

| File | What it does |
|---|---|
| `pssm.py` | `PSSM` class: for each sequence in a FASTA file, creates a temporary working folder, runs PSI-BLAST against a protein database (3 iterations, up to 5,000 alignments, E-value 0.001) and returns the ASCII PSSM file |
| `pssm_lib/blast.py` | `runPsiBlast()` wraps `psiblast`; `makeblastdb()` builds the database index if it is missing |
| `pssm_lib/workEnv.py` | `TemporaryEnv`: temporary folders and files for each job |
| `pssm_lib/datacache.py` | `DataCache`: caches PSSM files by the SHA-512 hash of the sequence, so a sequence is not searched twice |
| `pssm_lib/config.py` | `AAIDX10`: 10 normalised physicochemical property values for each amino acid |
| `pssm_class.ipynb`, `deneme.py` | Earlier drafts of the `PSSM` class |
| `ex.fasta`, `datasets/PP_matrix.fasta` | Example protein sequences |

## Usage

Requirements: Python 3, Biopython, NumPy, BLAST+ (`psiblast`, `makeblastdb`) and a protein FASTA database
(for example UniProt Swiss-Prot).

Set `fastafile` and `dbfile` at the top of `pssm.py`, then run:

```bash
python pssm.py
```

The PSSM files are written to `job.tmpd.*` folders in the current directory.

## Notes

- The `__pycache__/` folders are build files and can be deleted.
- `.vscode/launch.json` is an editor setting and is not needed to run the code.
