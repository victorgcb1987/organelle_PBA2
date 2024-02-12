# organelle_PBA2

Organelle PBA2 is a organelle genome assembler using Third Generation sequencing datasets like PacBIO or Oxford Nanopore. Following the original Organelle_PBA's (https://github.com/aubombarely/Organelle_PBA) premises it captures long reads from organelles by mapping reads to a reference genome.

The following software is required in order to run Organelle_PBA2:

# Requirements

Python3

Biopython (https://biopython.org/wiki/Download)

BLAST suite (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

FiltLong (https://github.com/rrwick/Filtlong)

Canu (https://github.com/marbl/canu)

Seqtk (https://github.com/lh3/seqtk)

Minimap2 (https://github.com/lh3/minimap2)

Racon (https://github.com/isovic/racon)

Organelle_PBA2 is designed to either check $PATH enviroment variable for the exectuables needed to run, but custom paths to executables can be defined like this

```
export MINIMAP2_PATH="/Path/to/executable/"
export NUCMER_PATH="/Path/to/executable/"
export FILTLONG_PATH="/Path/to/executable/"
export SEQTK_PATH="/Path/to/executable/"
export CANU_PATH="/Path/to/executable/"
export BLAST_PATH="/Path/to/executable/"
export RACON_PATH="/Path/to/executable/"
```

# Installation

Althought optional, It's strongly advised to install this program using a python virtualenv.

```
git clone https://github.com/victorgcb1987/organelle_PBA2.git
cd organelle_PBA2
python setup.py install
```

# Usage
This suite have three programs that can be used:

## Organelle assembly

```
python /organelle_PBA2/orgpba2/assemble_organelle.py --reference {reference_genome.fasta}  --sequence {reads.fq.gz} /
--threads {num_threads} --technology {pacbio|pacbio-hifi|nanopore}  --subsample_coverage {coverage} --out {output_dir} --curate_assembly
```
