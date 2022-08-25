#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python splitOneFastaIntoCertainNumber.py fastaFile outDir numberOfFiles
"""

import sys
import gzip

infileName = sys.argv[1]
outDir = sys.argv[2]
numberOfFiles = int(sys.argv[3])

def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)

numSeqsTotal = 0
with myopen(infileName, 'rb') as infile:
    for line in infile:
        if line[0] == ">":
            numSeqsTotal += 1

maxSeqPerFile = int(numSeqsTotal/float(numberOfFiles))+1

seqCounter = 0
outfile = myopen(outDir+'/split_'+str(seqCounter)+".fasta", "wb")
with myopen(infileName, 'rb') as infile:
    for line in infile:
        if line[0] == '>':
            seqCounter += 1
            if (seqCounter % maxSeqPerFile) == 0:
                outfile.close()
                outfile = myopen(outDir+'/split_'+str(seqCounter)+".fasta", "wb")
        outfile.write(line)
outfile.close()

