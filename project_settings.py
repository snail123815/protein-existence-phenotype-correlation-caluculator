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

REF = "ATMOS43"  # not used in step_0
EXP_DATA = Path("tested_strains_202403")  # Screening results

MBT_COLL_P = Path("MBT-collections/collective-faa/")
PROJ_PROTEOMES_P = Path("all_proteomes/")
DB_F = Path("all_proteomes_db") / "all_proteomes.fasta"
STRAINS_PICKLE = Path("step_0_gather_proteome_strains.pickle")
DOMTBLOUT_P = Path("domtblout.txt")
GATHER_DOMTBL_P = DOMTBLOUT_P.parent / (
    f"{DOMTBLOUT_P.stem}_E{str(GATHER_T_E)}"
    f"_DOME{str(GATHER_T_DOME)}_COV{str(GATHER_T_COV)}_LDIF{str(LEN_DIFF)}.tsv"
)
GATHER_MATCH_P = DOMTBLOUT_P.parent / (
    f"{DOMTBLOUT_P.stem}_matches_E{str(GATHER_T_E)}"
    f"_DOME{str(GATHER_T_DOME)}_COV{str(GATHER_T_COV)}_LDIF{str(LEN_DIFF)}.tsv"
)
PRESENCE_TABLE_P = Path(f"{GATHER_MATCH_P.stem}_presence.tsv")
