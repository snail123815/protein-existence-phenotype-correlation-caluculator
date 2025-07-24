import gzip
import pickle
from subprocess import run

from load_configs import (
    CONCATENATED_PROTEOMES_FILE,
    DOMTBLOUT_FILE,
    NCPU,
    STRAINS_PICKLE_FILE,
    T_DOME,
    T_E,
    T_INCDOME,
    T_INCE,
    TARGET_STRAIN,
)

REF_NEXT_LOC = 0
# Used if the program died in the middle. Use `get_location.py`
# to get where the program lefted out.

if __name__ == "__main__":
    with STRAINS_PICKLE_FILE.open("rb") as sp:
        phenotype_strains, all_strains = pickle.load(sp)

    for st in all_strains:
        if st == TARGET_STRAIN:
            ref_proteome_p = all_strains[st]
            break
    for phenotype, strains in phenotype_strains.items():
        for st in strains.items():
            pass

    with gzip.open(ref_proteome_p, "rt") as ref:
        ref.seek(REF_NEXT_LOC)
        jackhmmer_run = run(
            (
                f"jackhmmer -E {str(T_E)} --incE {str(T_INCE)} "
                f"--domE {str(T_DOME)}  --incdomE {str(T_INCDOME)} "
                f"--cpu {str(NCPU)} --domtblout {DOMTBLOUT_FILE} "
                f"- {CONCATENATED_PROTEOMES_FILE} 1>/dev/null"
            ),
            input=ref.read().encode(),
            shell=True,
            capture_output=True,
        )

    if jackhmmer_run.returncode != 0:
        print(jackhmmer_run.stderr.decode())
