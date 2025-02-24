from pathlib import Path

NCPU = 16

MIN_PROTEIN_LEN = 60
T_E = 1e-5
T_INCE = 1e-10
T_DOME = 1e-5
T_INCDOME = 1e-10
GATHER_T_E = 1e-10
GATHER_T_DOME = 1e-20
GATHER_T_COV = 0.7
LEN_DIFF = 0.2  # diff / min(qlen, tlen)

TARGET_STRAIN = "ATMOS43"
# Target strain, all of its proteins will be used, form rows in the
# presence table.

PHENOTYPE_TABLE_FILE = Path("tested_strains_202403")
# Contains the list of strains and their experimental data.
# Current format:
# MBT1|no conj.
# Target format:
# ID      Phenotype1      Phenotype2      ...
# MBT1    1    1    0

SOURCE_DATABASE_DIR = Path("MBT-collections/collective-faa/")
# Contains the target proteome databases. The tested strains in phenotype table
# should be present in this directory. File names are "*.faa.gz".

TEMP_PROTEOMICS_IN_TABLE_DIR = Path("all_proteomes/")
# Will generate this directory if not present. Symbolic links will be created
# pointing to the target proteome databases. As a result, it should contains
# the proteomes of all tested strains, including the target strain.

CONCATENATED_PROTEOMES_FILE = Path("all_proteomes_db") / "all_proteomes.fasta"
# Concatenated proteomes of all tested strains. Used as the database for
# jackhmmer. There should be no duplicate protein IDs in this file.

STRAINS_PICKLE_FILE = Path("step_0_gather_proteome_strains.pickle")
# Currently stores a dict of {phenotype: strains}

DOMTBLOUT_FILE = Path("domtblout.txt")
# Output of jackhmmer. Contains the domain hits of the target strain proteins.
# Will be used to generate the presence table.

# Some output files with generated file names based on the settings.
GATHER_DOMTBL_TSV = DOMTBLOUT_FILE.parent / (
    f"{DOMTBLOUT_FILE.stem}_E{str(GATHER_T_E)}"
    f"_DOME{str(GATHER_T_DOME)}_COV{str(GATHER_T_COV)}_LDIF{str(LEN_DIFF)}.tsv"
)
GATHER_MATCH_TSV = DOMTBLOUT_FILE.parent / (
    f"{DOMTBLOUT_FILE.stem}_matches_E{str(GATHER_T_E)}"
    f"_DOME{str(GATHER_T_DOME)}_COV{str(GATHER_T_COV)}_LDIF{str(LEN_DIFF)}.tsv"
)
PRESENCE_TSV = Path(f"{GATHER_MATCH_TSV.stem}_presence.tsv")
