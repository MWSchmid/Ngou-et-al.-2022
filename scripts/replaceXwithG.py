#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python replaceXwithG.py infileName outfileName mapForReversion
"""

import sys
import re
import gzip

try:
    infileName = sys.argv[1]
    outfileName = sys.argv[2]
    mapForReversion = sys.argv[3]
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

sequences = loadFasta(infileName)
replaceCollection = {}

numReplaced = 0
with myopen(outfileName, 'w') as outfile:
    for name, seq in sequences.items():
        seq = seq.upper()
        seqList = list(seq)
        for i, c in enumerate(seq):
            if c in ["X", "B", "J", "Z", "-", "U", ".", "*"]:
                numReplaced += 1
                seqList[i] = "G"
                entry = '\t'.join([name, str(i), c])
                try:
                    replaceCollection[name].append(entry)
                except KeyError:
                    replaceCollection[name] = [entry]
        outSeq = ''.join(seqList)
        outfile.write(''.join([">", name, "\n", outSeq, "\n"]))

with myopen(mapForReversion, 'w') as outfile:
    for name, changes in replaceCollection.items():
        outfile.write('\n'.join(changes)+'\n')

print >> sys.stderr, "replaced %d entries" % numReplaced
