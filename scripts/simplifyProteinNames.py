#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python simplifyProteinNames.py prefix infile outfile mapfile
"""

import sys
import gzip

try:
    prefix = sys.argv[1]
    fastaInfile = sys.argv[2]
    fastaOutfile = sys.argv[3]
    mapOutfile = sys.argv[4]
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

with myopen(mapOutfile, "w") as mapOut, myopen(fastaOutfile, "w") as fastaOut:
    seqCounter = 0
    for name, seq in sequences.items():
        seqCounter += 1
        newName = ''.join([prefix, str(seqCounter)])
        outString = ''.join(['>', newName, '\n', seq, '\n'])
        fastaOut.write(outString)
        outString = ''.join([newName, "|||", name, '\n'])
        mapOut.write(outString)
        
        
        

