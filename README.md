# TCRseq
A tool for building the TCR sequence, given V-gene name, J-gene name and CDR3 sequence.

The data directory contains the fasta files, holding amino acid sequences of the V and J genes for alpha(A) and beta(B) TCR subunits.
These have been downloaded from IMGT (https://www.imgt.org/vquest/refseqh.html) and should be updated whenever this database gets an update.

## dependencies
 * biopython==1.83

## usage:
From this repository, run:

```
python3 tcrseq.py VGENE_NAME JGENE_NAME CDR3_SEQUENCE
```

This returns the full TCR sequence.

This python file can also be used as a module:

```
from tcrseq import get_full_tcr_sequence

full_sequence = get_full_sequence("TRBV1", "TRBJ1", "CTSSQAEFAFANTEA")
```
