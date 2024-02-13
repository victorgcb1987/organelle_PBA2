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
--threads {num_threads} --technology {pacbio|pacbio-hifi|nanopore}  --subsample_coverage {coverage} /
--out {output_dir} --curate_assembly
```
Plastid/Mitchondria reference genome (--reference), long reads dataset (--sequence) and output dir (--out) are the only arguments required to run the assembly.
Command above shows the typical calling of this program but there is other optional arguments that can be declared, like specifiying how minimap2, filtlong or 
canu options (--minimap2, --canu_options, -filtlong_options).

The --curate_assembly will try to curate the assembly by homology against the refence genome provided in the command call, so it's advisable to use a reference genome from the same genus.

Assemble organelle stores each step required in order to run the assembly in it's own folder, so following rerunnings of this programs skips already done steps.
You can force to redo each step declaring the following, optional arguments when rerunning the program:
```
--force_mapping
--force_extract_mapped_reads
--force_assembly
--force_subsampling
```

## Calculate heteroplasmy

```
python calculate_heteroplasmy.py --assembly {reference_genome.fasta} --sequences {reads.fq.gz} --out {output_dir} /
--threads {num_threads}  --technology {pacbio|pacbio-hifi|nanopore}
```


## NUPTs detection

```
python indentify_nuclear_insertions.py --nuclear_assembly {nuclear_reference.fasta} \
--organelle_assembly {organelle_assembly.fasta} --exclude {organelle_to_compare.fasta} \
--technology {pacbio|pacbio-hifi|nanopore}  --threads {num_threads} \
--sequences {reads.fq.gz} --length {int} --out {out_dir}
```
