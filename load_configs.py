from pathlib import Path

import yaml

with open(Path(__file__).parent / "config_jackhmmer_thresholds.yaml") as f:
    thresholds_config = yaml.safe_load(f)

with open(Path(__file__).parent / "config_project.yaml") as f:
    project_config = yaml.safe_load(f)

MIN_PROTEIN_LEN = int(thresholds_config["MIN_PROTEIN_LEN"])
T_E = float(thresholds_config["T_E"])
T_INCE = float(thresholds_config["T_INCE"])
T_DOME = float(thresholds_config["T_DOME"])
T_INCDOME = float(thresholds_config["T_INCDOME"])
GATHER_T_E = float(thresholds_config["GATHER_T_E"])
GATHER_T_DOME = float(thresholds_config["GATHER_T_DOME"])
GATHER_T_COV = float(thresholds_config["GATHER_T_COV"])
LEN_DIFF = float(thresholds_config["LEN_DIFF"])

NCPU = int(project_config["NCPU"])
TARGET_STRAIN = str(project_config["TARGET_STRAIN"])
PHENOTYPE_TABLE_FILE = Path(project_config["PHENOTYPE_TABLE_FILE"])
SOURCE_DATABASE_DIR = Path(project_config["SOURCE_DATABASE_DIR"])
TEMP_PROTEOMICS_IN_TABLE_DIR = Path(project_config["TEMP_PROTEOMICS_IN_TABLE_DIR"])
CONCATENATED_PROTEOMES_FILE = Path(project_config["CONCATENATED_PROTEOMES_FILE"])
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
PRESENCE_TSV = Path(GATHER_MATCH_TSV.parent / f"presence_{GATHER_MATCH_TSV.stem}.tsv")
