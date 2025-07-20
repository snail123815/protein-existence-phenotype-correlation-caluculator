# Read tested_strains TSV file with headers, get phenotype lists dynamically:
# Format: ID\tPhenotype1\tPhenotype2\t...
#         strain1\t1\t0\t...
#         strain2\t1\t1\t...
# Then, make a database containing only *all tested strains*

import gzip
import pickle
import re
from pathlib import Path

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


def handle_missing_proteome(strain_name):
    """
    Handle missing proteome files by asking user for input.

    Args:
        strain_name (str): Name of the strain with missing proteome

    Returns:
        bool: True to continue, False to exit
    """
    print(f"\nWARNING: Proteome of strain '{strain_name}' not found.")
    while True:
        choice = (
            input("Do you want to continue without this strain? (y/n/a): ")
            .lower()
            .strip()
        )
        if choice in ["y", "yes"]:
            return True
        elif choice in ["n", "no"]:
            print("Exiting program...")
            return False
        elif choice in ["a", "all"]:
            print("Continuing with all missing strains...")
            return "all"
        else:
            print("Please enter 'y' (yes), 'n' (no), or 'a' (all missing).")


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
    df = pd.read_csv(
        PHENOTYPE_TABLE_FILE,
        sep="\t",
        index_col=0,
    )
    df = df.astype(float)

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

    print(f"Total strains in {PHENOTYPE_TABLE_FILE}: {df.shape[0]}")
    print(f"Phenotypes found: {phenotype_names}")
    for phenotype, strains in phenotype_strains.items():
        print(f"  {phenotype}: {len(strains)} strains")

    # Find proteome files for each strain in each phenotype
    continue_all_missing = False  # Flag to auto-continue for all missing files
    all_strains: dict[str, Path] = {}  # Store all strains for later use
    for st in list(df.index):
        found = False
        for f in SOURCE_DATABASE_DIR.glob("*.faa.gz"):
            if st in re.split(r"_|\.", f.name):
                # print(f'Strain {st} path {f}')
                found = True
                # Create a symlink in TEMP_PROTEOMICS_IN_TABLE_DIR
                symlink_path = TEMP_PROTEOMICS_IN_TABLE_DIR / f.name
                try:
                    symlink_path.symlink_to(
                        f.relative_to(
                            TEMP_PROTEOMICS_IN_TABLE_DIR, walk_up=True
                        )
                    )
                except FileExistsError:
                    pass
                all_strains[st] = symlink_path
                break
        if not found:
            continue_anyway = True
            if not continue_all_missing:
                continue_anyway = handle_missing_proteome(st)
                if continue_anyway is False:
                    # User chose to exit
                    exit(1)
                elif continue_anyway == "all":
                    continue_all_missing = True
            if continue_anyway or continue_all_missing:
                print(f"Skipping strain {st} (missing proteome).")
                all_strains.pop(st, None)
                for phenotype, strains in phenotype_strains.items():
                    if st in strains:
                        phenotype_strains[phenotype].pop(st, None)

    CONCATENATED_PROTEOMES_FILE.parent.mkdir(exist_ok=True)

    print("Making database fasta:")
    if CONCATENATED_PROTEOMES_FILE.exists():
        print(f"Removing existing {CONCATENATED_PROTEOMES_FILE}.")
        CONCATENATED_PROTEOMES_FILE.unlink()
    CONCATENATED_PROTEOMES_FILE.touch()
    with CONCATENATED_PROTEOMES_FILE.open("at", encoding="utf-8") as db_handle:
        for st, proteome_p in tqdm(all_strains.items()):
            with gzip.open(proteome_p, "rt") as source:
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
