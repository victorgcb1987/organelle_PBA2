

from importlib.util import MAGIC_NUMBER

EXECUTABLES_REQUIREMENTS = {"minimap2": {"executable": "minimap2",
                                         "user_path": "MINIMAP2_PATH"},
                            "blastn": {"executable": "blastn",
                                       "user_path": "BLAST_PATH"},
                            "makeblastdb": {"executable": "makeblastdb",
                                            "user_path": "BLAST_PATH"},
                            "seqtk": {"executable": "seqtk",
                                      "user_path": "SEQTK_PATH"},
                            "canu": {"executable": "canu",
                                     "user_path": "CANU_PATH"},
                            "filtlong": {"executable": "filtlong",
                                     "user_path": "FILTLONG_PATH"},
                            "racon": {"executable": "racon",
                                      "user_path": "RACON_PATH"},
                            "nucmer": {"executable": "nucmer",
                                      "user_path": "NUCMER_PATH"}}

MAGIC_NUMS_COMPRESSED = [b'\x1f\x8b\x08', b'\x42\x5a\x68', 
                         b'\x50\x4b\x03\x04']

OUTPUT_FOLDERS = {"haplotypes": "07_haplotypes",
                  "minimap2": "01_mapping_minimap2",
                  "mapped_seqs": "02_mapped_seqs_seqtk_filtlong",
                  "subsampled_seqs": "02_mapped_seqs_seqtk_filtlong",
                  "canu": "03_assembly_canu",
                  "blast_db": "00_blasts",
                  "polished_assembly": "04_polished_assembly_racon",
                  "curated_assembly": "06_curated_assembly_mummer",
                  "no_redundant": "05_remove_circularity_redundancy",
                  "repetitive": "05_rescaffold_repetitive_regions",
                  }
OUTPUT_FILENAMES = {"haplotypes": "06_haplotypes.fasta",
                    "minimap2": "01_mapping_minimap2.paf",
                    "mapped_seqs": "02_mapped_seqs.fq.gz",
                    "subsampled_seqs": "02_subsampled_seqs_cov{}.fq.gz",
                    "assembly": "assembly.fasta",
                    "curated_assembly": "curated_assembly.fasta",
                    "curated_contigs": "curated_contigs.fasta",
                    "contigs_without_overlaps": "contigs_without_overlaps.fasta",
                    "concatenated_contigs": "concatenated_contigs.fasta",
                    "polished_assembly": "04_polished_assembly_iteration_{}.fasta",
                    "mapping_for_polishing": "04_polishing_mapping_iteration_{}.paf",
                    "mapping_for_haplotypes": "07_haplotypes_reads_mapped.paf",
                    "oriented_contigs": "oriented_contigs.fasta",
                    "no_redundant": "05_non_redundant.fasta"}

COMPLEMENTARY_NUCLEOTIDE = {"A": "T", "T": "A", 
                            "C": "G", "G": "C", 
                            "N": "N"}

