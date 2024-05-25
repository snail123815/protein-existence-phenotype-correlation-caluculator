# If jackhmmer died from memory failer, use this to query the location of the
# input file where it should start in the middle.

import gzip
from pathlib import Path


target = ">TMLOC_02642"
ref_proteome_p = Path(
    "./MBT-collection/collective-faa/Streptomyces_sp._ATMOS43.faa.gz"
)

# Check the first appearance
with gzip.open(ref_proteome_p, "rt") as rpp:
    assert rpp.seekable()

    l = rpp.readline()
    while l:
        if l.startswith(target):
            print(rpp.tell())
            break
        l = rpp.readline()

# Find the next gene location
with gzip.open(ref_proteome_p, "rt") as rpp:
    rpp.seek(964405)
    l = rpp.readline()
    old_loc = 0
    while l:
        print(l)
        current_loc = rpp.tell()
        if l.startswith(">"):
            print(old_loc)
            break
        old_loc = current_loc
        l = rpp.readline()

# Check the result
with gzip.open(ref_proteome_p, "rt") as rpp:
    rpp.seek(964563)
    for i in range(5):
        print(rpp.readline())
