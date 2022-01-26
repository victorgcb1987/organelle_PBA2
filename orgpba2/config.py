

from importlib.util import MAGIC_NUMBER


EXECUTABLES_REQUIREMENTS = {"minimap2": {"executables": ["minimap2"],
                                     "user_path": "MINIMAP2_PATH"},
                         "blast": {"executables": ["blastn", "makeblastdb"],
                                   "user_path": "BLAST_PATH"},
                                   }

MINIMAP2_LONGARGS = {"indexer": ["idx-no-seq", "alt", 
                                 "alt-drop"],
                     "mapper": ["q-occ-frac", "dual",
                                "rmq", "hard-mask-level",
                                "mask-len", "max-chain-skip",
                                "max-chain-iter", "chain-gap-scale",
                                "no-long-join", "splice", "sr",
                                "split-prefix", "frag", "for-only",
                                "rev-only", "heap-sort", "no-pairing"],
                    "align": ["end-bonus", "score-N", "splice-flank",
                              "junc-bed", "junc-bonus", "end-seed-pen",
                              "no-end-flt", "cap-sw-mem", "cap-kalloc"],
                    "IO": ["cs", "MD", "eqx", "seed", "secondary",
                           "max-qlen", "paf-no-hit", "sam-hit-only",
                            "version"],
                    "misc": ["no-kalloc", "print-qname", "print-seeds"]}

MAGIC_NUMS_COMPRESSED = [b'\x1f\x8b\x08', b'\x42\x5a\x68', 
                         b'\x50\x4b\x03\x04']

OUTPUT_FOLDERS = {"minimap2": "01_mapping_minimap2", 
                  "canu": "02_canu"}