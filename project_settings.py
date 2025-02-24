import json
from pathlib import Path

with open(Path(__file__).parent / "config_jackhmmer_thresholds.json") as f:
    thresholds_config = json.load(f)

with open(Path(__file__).parent / "config_project.json") as f:
    project_config = json.load(f)

MIN_PROTEIN_LEN = thresholds_config["MIN_PROTEIN_LEN"]
T_E = thresholds_config["T_E"]
T_INCE = thresholds_config["T_INCE"]
T_DOME = thresholds_config["T_DOME"]
T_INCDOME = thresholds_config["T_INCDOME"]
GATHER_T_E = thresholds_config["GATHER_T_E"]
GATHER_T_DOME = thresholds_config["GATHER_T_DOME"]
GATHER_T_COV = thresholds_config["GATHER_T_COV"]
LEN_DIFF = thresholds_config["LEN_DIFF"]

NCPU = project_config["NCPU"]
TARGET_STRAIN = project_config["TARGET_STRAIN"]
PHENOTYPE_TABLE_FILE = Path(project_config["PHENOTYPE_TABLE_FILE"])
SOURCE_DATABASE_DIR = Path(project_config["SOURCE_DATABASE_DIR"])
TEMP_PROTEOMICS_IN_TABLE_DIR = Path(
    project_config["TEMP_PROTEOMICS_IN_TABLE_DIR"]
)
CONCATENATED_PROTEOMES_FILE = Path(
    project_config["CONCATENATED_PROTEOMES_FILE"]
)
STRAINS_PICKLE_FILE = Path(project_config["STRAINS_PICKLE_FILE"])
DOMTBLOUT_FILE = Path(project_config["DOMTBLOUT_FILE"])

GATHER_DOMTBL_TSV = DOMTBLOUT_FILE.parent / (
    f"{DOMTBLOUT_FILE.stem}_E{str(GATHER_T_E)}"
    f"_DOME{str(GATHER_T_DOME)}_COV{str(GATHER_T_COV)}_LDIF{str(LEN_DIFF)}.tsv"
)
GATHER_MATCH_TSV = DOMTBLOUT_FILE.parent / (
    f"{DOMTBLOUT_FILE.stem}_matches_E{str(GATHER_T_E)}"
    f"_DOME{str(GATHER_T_DOME)}_COV{str(GATHER_T_COV)}_LDIF{str(LEN_DIFF)}.tsv"
)
PRESENCE_TSV = Path(f"{GATHER_MATCH_TSV.stem}_presence.tsv")
