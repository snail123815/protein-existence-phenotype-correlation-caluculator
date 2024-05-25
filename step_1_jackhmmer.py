import pickle
import gzip
from pathlib import Path
from subprocess import run

from step_0_gather_proteome import DB_F, STRAINS_PICKLE, REF

T_E = 1e-5
T_INCE = 1e-10
T_DOME = 1e-5
T_INCDOME = 1e-10
NCPU = 16

DOMTBLOUT_P = Path("domtblout.txt")
REF_NEXT_LOC = 0
# Used if the program died in the middle. Use `get_location.py`
# to get where the program lefted out.

if __name__ == "__main__":
    with STRAINS_PICKLE.open('rb') as sp:
        strains_paths_dict = pickle.load(sp)

    strains_doubleconj = strains_paths_dict["Double conj."]
    strains_singleconj = strains_paths_dict["Single conj."]
    strains_noconj = strains_paths_dict["No conj."]

    ref_proteome_p = strains_doubleconj[REF]
    with gzip.open(ref_proteome_p, 'rt') as ref:
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
