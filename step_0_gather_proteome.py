# Read tested_strains TSV file with headers, get phenotype lists dynamically:
# Format: ID\tPhenotype1\tPhenotype2\t...
#         strain1\t1\t0\t...
#         strain2\t1\t1\t...
# Then, make a database containing only *all tested strains*

import gzip
import pickle
import re

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from load_configs import (
    CONCATENATED_PROTEOMES_FILE,
    MIN_PROTEIN_LEN,
    PHENOTYPE_TABLE_FILE,
    SOURCE_DATABASE_DIR,
    STRAINS_PICKLE_FILE,
    TEMP_PROTEOMICS_IN_TABLE_DIR,
)

# STRAINS_PICKLE = Path("step_0_gather_proteome_strains.pickle")
# Saves a dict of dict:
# {
#    "Phenotype1": {"strain1": PosixPath("..."), "strain2": PosixPath("...")},
#    "Phenotype2": {"strain3": PosixPath("..."), "strain4": PosixPath("...")},
#    ...
# }
# With each strain mapped to its proteome file path

if __name__ == "__main__":
    # Dictionary to store strains for each phenotype
    # Structure: {"Phenotype1": {"strain1": "", "strain2": ""}, ...}
    phenotype_strains = {}

    # Remove all files in TEMP_PROTEOMICS_IN_TABLE_DIR
    if TEMP_PROTEOMICS_IN_TABLE_DIR.exists():
        for f in TEMP_PROTEOMICS_IN_TABLE_DIR.iterdir():
            f.unlink()
    else:
        TEMP_PROTEOMICS_IN_TABLE_DIR.mkdir()

    # Read the TSV file using pandas
    df = pd.read_csv(PHENOTYPE_TABLE_FILE, sep="\t", index_col=0, dtype=float)

    # Get phenotype names from column headers
    phenotype_names = df.columns.tolist()

    # Initialize dictionaries for each phenotype
    for phenotype in phenotype_names:
        phenotype_strains[phenotype] = {}

    # Process each phenotype column
    for phenotype in phenotype_names:
        # Get strains where the phenotype value is > 0 (and not NaN)
        positive_strains = df[df[phenotype] > 0][phenotype]
        for strain, value in positive_strains.items():
            phenotype_strains[phenotype][strain] = value

    print(f"Total strains in {PHENOTYPE_TABLE_FILE}: {len(df)}")
    print(f"Phenotypes found: {phenotype_names}")
    for phenotype, strains in phenotype_strains.items():
        print(f"  {phenotype}: {len(strains)} strains")

    # Find proteome files for each strain in each phenotype
    for phenotype_name, strains in phenotype_strains.items():
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
    for phenotype_name, strains in phenotype_strains.items():
        with CONCATENATED_PROTEOMES_FILE.open("at") as db_handle:
            for st in tqdm(strains, desc=phenotype_name):
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
        pickle.dump(phenotype_strains, sp)
