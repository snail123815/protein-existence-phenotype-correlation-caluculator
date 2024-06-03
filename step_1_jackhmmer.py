import gzip
import pickle
from subprocess import run

from project_settings import (DB_F, DOMTBLOUT_P, NCPU, REF, STRAINS_PICKLE,
                              T_DOME, T_E, T_INCDOME, T_INCE)

REF_NEXT_LOC = 0
# Used if the program died in the middle. Use `get_location.py`
# to get where the program lefted out.

if __name__ == "__main__":
    with STRAINS_PICKLE.open("rb") as sp:
        strains_paths_dict = pickle.load(sp)

    strains_doubleconj = strains_paths_dict["Double conj."]
    strains_singleconj = strains_paths_dict["Single conj."]
    strains_noconj = strains_paths_dict["No conj."]

    ref_proteome_p = strains_doubleconj[REF]
    with gzip.open(ref_proteome_p, "rt") as ref:
        ref.seek(REF_NEXT_LOC)
        jackhmmer_run = run(
            (
                f"jackhmmer -E {str(T_E)} --incE {str(T_INCE)} "
                f"--domE {str(T_DOME)}  --incdomE {str(T_INCDOME)} "
                f"--cpu {str(NCPU)} --domtblout {DOMTBLOUT_P} "
                f"- {DB_F} 1>/dev/null"
            ),
            input=ref.read().encode(),
            shell=True,
            capture_output=True,
        )

    if jackhmmer_run.returncode != 0:
        print(jackhmmer_run.stderr.decode())
