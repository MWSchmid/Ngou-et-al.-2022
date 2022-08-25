#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python revertSequences.py infileName mapForReversion
"""

import sys
import re
import gzip

try:
    infileName = sys.argv[1]
    mapForReversion = sys.argv[2]
except IndexError:
    print >> sys.stderr, __doc__
    sys.exit(1)

def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)

replaceCollection = {}
with myopen(mapForReversion) as infile:
    for line in infile:
        name, pos, base = line.strip().split('\t')
        pos = int(pos)
        try:
            replaceCollection[name].append([pos, base])
        except KeyError:
            replaceCollection[name] = [[pos, base]]

numReplaced = 0
with myopen(infileName) as infile:
    for line in infile:
        name, start, end, seq = line.strip().split('\t')
        if name not in replaceCollection:
            print line.strip()
            continue
        start = int(start)-1
        end = int(end)-1
        if end < start:
            print >> sys.stderr, "weird"
        seqList = list(seq)
        for pos, base in replaceCollection[name]:
            if pos > start and pos <= end:
                posInSeq = pos - start - 1  # that was trial and error
                if seqList[posInSeq] != "G":
                    print >> sys.stderr, pos
                    print >> sys.stderr, start
                    print >> sys.stderr, end
                    print >> sys.stderr, posInSeq
                    print >> sys.stderr, seqList
                else:
                    numReplaced += 1
                    seqList[posInSeq] = base
            else:
                continue
        seq = ''.join(seqList)
        print '\t'.join([name, str(start+1), str(end+1), seq])
        
print >> sys.stderr, "replaced %d entries" % numReplaced
