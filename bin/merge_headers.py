#!/usr/bin/env python3

import os
import re
import subprocess
import sys
from signal import signal, SIGPIPE, SIG_DFL

"""
where SN field in SQ header record matches one in the given dict file,
replace that SQ record with that in the dict file,
but propagate any 'AH' fields provided by the aligner"
"""


def exitWithError(error):
    print("Error: {}".format(error))
    sys.exit(2)


# -------------------------------------------------


def querySamtools(target):
    result = subprocess.run(
        ["samtools", "view", "-H", target],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode == 0:
        return result.stdout.decode("UTF-8")
    else:
        print("Error: {}".format(result.stderr.decode("UTF-8")))
        return None


# -------------------------------------------------


def getHeaderValue(line, target):
    if target in line:
        start = line.index(target)
        end = line.index("\t", start)
        return line[start:end]
    return ""


# -------------------------------------------------


def main():
    # Program arguments:
    # 1. Original cram file with adapters clipped
    # 2. New mapped bam
    # 3. Reference file path
    if len(sys.argv) < 3:
        exitWithError("Missing mandatory arguments")

    signal(SIGPIPE, SIG_DFL)
    regex = re.compile("@SQ.*\tSN:([^\t]+)")

    old = querySamtools(sys.argv[1])
    new = querySamtools(sys.argv[2])

    lastPGID = ""
    filteredOldArray = []
    for l in old.split("\n"):

        filteredOldArray.append(l)
        if l.startswith("@PG") and "ID:" in l:
            lastPGID = getHeaderValue(l, "ID:").replace("ID", "\tPP")

    old = "\n".join(filteredOldArray)

    lookup = {}
    # Dictionary file should be named like so: 'path/reference.dict'
    dictionary_file = os.path.splitext(sys.argv[3])[0]+".dict"
    with open(dictionary_file, "r") as fd:
        for line in fd:
            if "\t" in line:
                refBits = line.split("\t")
                lookup[refBits[1]] = line.strip()

    filteredNewArray = []
    newPGs = []
    newlines = new.split("\n")
    for line in newlines:
        if regex.search(line):
            refBits = line.split("\t")
            if refBits[1] in lookup:
                filteredNewArray.append(lookup[refBits[1]])
        else:
            if "@PG" in line:
                if "PP:" not in line:
                    newPGs.append(line + lastPGID)
                else:
                    newPGs.append(line)
            else:
                filteredNewArray.append(line)

    new = "\n".join(filteredNewArray)

    sys.stdout.write(old.strip())
    sys.stdout.write("\n")
    sys.stdout.write("\n".join(reversed(newPGs)).strip())
    sys.stdout.write("\n")
    sys.stdout.write(new.strip())
    sys.stdout.write("\n")


# -------------------------------------------------
if __name__ == "__main__":
    main()
else:
    print(__name__)


