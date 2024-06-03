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

from project_settings import (DB_F, EXP_DATA, MBT_COLL_P, MIN_PROTEIN_LEN,
                              PROJ_PROTEOMES_P, STRAINS_PICKLE)

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

    for f in PROJ_PROTEOMES_P.iterdir():
        f.unlink()

    with EXP_DATA.open() as expdata:
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

        print(f"Total lines in {EXP_DATA}: {i+1}")

    for strains in [strains_doubleconj, strains_noconj, strains_singleconj]:
        for st in list(strains.keys()):
            found = False
            for f in MBT_COLL_P.glob("*.faa.gz"):
                if st in re.split(r"_|\.", f.name):
                    strains[st] = f
                    # print(f'Strain {st} path {f}')
                    found = True
                    (PROJ_PROTEOMES_P / f.name).symlink_to(
                        f.relative_to(PROJ_PROTEOMES_P, walk_up=True)
                    )
                    break
            if not found:
                print(f"Proteome of strain {st} not found.")
                strains.pop(st)

    DB_F.parent.mkdir(exist_ok=True)

    print("Making database fasta:")
    for phenotype, strains in zip(
        ["Double conj.", "Single conj.", "No conj."],
        [strains_doubleconj, strains_singleconj, strains_noconj],
    ):
        with DB_F.open("at") as db_handle:
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

    print(f"Database fasta file {DB_F}.")
    if STRAINS_PICKLE.exists():
        STRAINS_PICKLE.unlink()
    STRAINS_PICKLE.touch()
    with STRAINS_PICKLE.open("wb") as sp:
        pickle.dump(
            {
                "Double conj.": strains_doubleconj,
                "Single conj.": strains_singleconj,
                "No conj.": strains_noconj,
            },
            sp,
        )
