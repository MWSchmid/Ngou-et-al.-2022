#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python convertPhyloXmlToNewick.py infileName outfileName
"""

import sys
from Bio import Phylo

try:
    infileName = sys.argv[1]
    outfileName = sys.argv[2]
except IndexError:
    print >> sys.stderr, __doc__
    sys.exit(1)

Phylo.convert(infileName, "phyloxml", outfileName, "newick")
