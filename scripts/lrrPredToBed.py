#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python lrrPredToBed.py infileName
"""

import sys
import re
import gzip

try:
    infileName = sys.argv[1]
except IndexError:
    print >> sys.stderr, __doc__
    sys.exit(1)

def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)

numWarnings = 0
numEmpty = 0
curName = "ukn"
with myopen(infileName, 'r') as infile:
    for line in infile:
        fields = line.strip().split()
        try:
            if fields[0] == "Warning:":
                numWarnings += 1
                continue
            if fields[0] == "Prediction":
                curName = fields[4].replace(':', '')
                continue
        except IndexError:
            numEmpty += 1
            continue
        seq = fields[3].replace(',', '')
        start = int(fields[2].replace(',', ''))
        end = start + len(seq)
        print '\t'.join([curName, str(start), str(end), seq])

print >> sys.stderr, "had %d warnings and %d empty lines" % (numWarnings, numEmpty)
