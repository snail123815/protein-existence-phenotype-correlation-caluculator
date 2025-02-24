# Read tested_strains file, get 3 lists:
# 1. Double conj.
# 2. Single conj.
# 3. No conj.
# Then, make a database containing only *all tested strains*

import gzip
import pickle
import re
from pathlib import Path

from Bio import SeqIO
from tqdm import tqdm

from project_settings import (
    CONCATENATED_PROTEOMES_FILE,
    PHENOTYPE_TABLE_FILE,
    SOURCE_DATABASE_DIR,
    MIN_PROTEIN_LEN,
    TEMP_PROTEOMICS_IN_TABLE_DIR,
    STRAINS_PICKLE_FILE,
)

# STRAINS_PICKLE = Path("step_0_gather_proteome_strains.pickle")
# Saves a dict of dict:
# {
#    "Double conj.": strains_doubleconj,
#    "Single conj.": strains_singleconj,
#    "No conj.": strains_noconj,
# }
# With each {"MBT1": PosixPath("MBT-collection/collective-faa/MBT1.fa.gz")}

if __name__ == "__main__":
    strains_noconj = {}
    # {"MBT1": PosixPath("MBT-collection/collective-faa/MBT1.fa.gz")}
    strains_singleconj = {}
    strains_doubleconj = {}

    # Remove all files in TEMP_PROTEOMICS_IN_TABLE_DIR
    if TEMP_PROTEOMICS_IN_TABLE_DIR.exists():
        for f in TEMP_PROTEOMICS_IN_TABLE_DIR.iterdir():
            f.unlink()
    else:
        TEMP_PROTEOMICS_IN_TABLE_DIR.mkdir()

    with PHENOTYPE_TABLE_FILE.open() as expdata:
        for i, l in enumerate(expdata):
            strain, exp = l.split("|")
            strain = strain.strip()
            exp = exp.strip()
            match exp:
                case "no conj.":
                    strains_noconj[strain] = ""
                case "Single AA conj.":
                    strains_singleconj[strain] = ""
                case "di-peptide conj.":
                    strains_doubleconj[strain] = ""
                case _:
                    print(f"Case {exp} not known")

        print(f"Total lines in {PHENOTYPE_TABLE_FILE}: {i+1}")

    for strains in [strains_doubleconj, strains_noconj, strains_singleconj]:
        for st in list(strains.keys()):
            found = False
            for f in SOURCE_DATABASE_DIR.glob("*.faa.gz"):
                if st in re.split(r"_|\.", f.name):
                    strains[st] = f
                    # print(f'Strain {st} path {f}')
                    found = True
                    (TEMP_PROTEOMICS_IN_TABLE_DIR / f.name).symlink_to(
                        f.relative_to(
                            TEMP_PROTEOMICS_IN_TABLE_DIR, walk_up=True
                        )
                    )
                    break
            if not found:
                print(f"Proteome of strain {st} not found.")
                strains.pop(st)

    CONCATENATED_PROTEOMES_FILE.parent.mkdir(exist_ok=True)

    print("Making database fasta:")
    for phenotype, strains in zip(
        ["Double conj.", "Single conj.", "No conj."],
        [strains_doubleconj, strains_singleconj, strains_noconj],
    ):
        with CONCATENATED_PROTEOMES_FILE.open("at") as db_handle:
            for st in tqdm(strains, desc=phenotype):
                proteome_p = strains[st]
                with gzip.open(strains[st], "rt") as source:
                    prots = []
                    for prot in SeqIO.parse(source, "fasta"):
                        prot.id = f"{st}_{prot.id}"
                        if len(prot) < MIN_PROTEIN_LEN:
                            continue
                        prots.append(prot)
                    SeqIO.write(prots, db_handle, "fasta")

    print(f"Database fasta file {CONCATENATED_PROTEOMES_FILE}.")
    if STRAINS_PICKLE_FILE.exists():
        STRAINS_PICKLE_FILE.unlink()
    STRAINS_PICKLE_FILE.touch()
    with STRAINS_PICKLE_FILE.open("wb") as sp:
        pickle.dump(
            {
                "Double conj.": strains_doubleconj,
                "Single conj.": strains_singleconj,
                "No conj.": strains_noconj,
            },
            sp,
        )
