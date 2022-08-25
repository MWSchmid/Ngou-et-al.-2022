#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python extractTMproteins.py.py maxTM tmhmmRed fastaFileIn fastaFileOut
"""

import sys
import re
import gzip

try:
    maxTM = int(sys.argv[1])
    tmhmmFile = sys.argv[2]
    fastaInfile = sys.argv[3]
    fastaOutfile = sys.argv[4]
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

sequences = loadFasta(fastaInfile)

hasTMMs = set([])
with myopen(tmhmmFile) as infile:
    for line in infile:
        line = line.strip()
        if line[0] != "#":
            continue
        fields = line.split(' ')
        if fields[2] != "Number":
            continue
        if int(fields[-1]) > maxTM:
            name = fields[1]
            hasTMMs.add(name)

counter = 0
with myopen(fastaOutfile, 'w') as outfile:
    for name, seq in sequences.items():
        if name in hasTMMs:
            counter += 1
            outString = ''.join(['>', name, '\n', seq, '\n'])
            outfile.write(outString)

print >> sys.stderr, "found %d peptides." % (counter)

