#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python searchLongerProteins.py minLen fastaFile
"""

import sys
import re
import gzip

try:
    minLen = int(sys.argv[1])
    fastaFile = sys.argv[2]
except IndexError:
    print >> sys.stderr, __doc__
    sys.exit(1)

def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)

def loadFasta(fastaFile):
    genomeDict = {}
    with myopen(fastaFile, 'r') as infile:
        for line in infile:
            if line[0] == '>':
                chromosome = line.strip().split(' ')[0][1:]
                genomeDict.update({chromosome:[]})
            else:
                genomeDict[chromosome].append(line[0:-1])
        for chromosome in genomeDict:
            sequence = ''.join(genomeDict[chromosome])
            genomeDict[chromosome] = sequence
    return genomeDict

sequences = loadFasta(fastaFile)
counter = 0
for name, seq in sequences.items():
    if len(seq) >= minLen:
        print '>'+name
        print seq
        counter += 1


print >> sys.stderr, "found %d peptides." % (counter)
