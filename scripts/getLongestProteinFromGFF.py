# -*- coding: utf-8 -*-
"""
python getLongestProteinFromGFF.py proteinFasta gffInfileName outFasta
"""

import sys
import re
import gzip

try:
    proteinFasta = sys.argv[1]
    gffInfileName = sys.argv[2]
    outfileName = sys.argv[3]
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
        useTranscript = False
        for line in infile:
            if line[0] == '>':
                fields = line.strip()[1:].split()
                if len(fields) == 0:
                    chromosome = fields[0]
                else:
                    if useTranscript:
                        chromosome = fields[transcriptField]
                        chromosome = re.sub("transcript=", "", chromosome)
                    else:
                        for i,curField in enumerate(fields):
                            if "transcript=" in curField:
                                useTranscript = True
                                transcriptField = i
                        if useTranscript:
                            chromosome = fields[transcriptField]
                            chromosome = re.sub("transcript=", "", chromosome)
                        else:
                            chromosome = fields[0]
                # some have suffixes, hope that the numbers match...
                if chromosome[:2] == "Am":
                    chromosome = re.sub("\.P", ".T", chromosome)
                if chromosome[:5] == "AALBA":
                    chromosome = re.sub("P", "T", chromosome)
                #if chromosome[:4] == "VITM":
                #    chromosome = re.sub("\.p", ".T", chromosome)
                chromosome = re.sub("\-P$", "", chromosome)
                chromosome = re.sub("\.p$", "", chromosome)
                chromosome = re.sub("\.cds\d$", "", chromosome)
                chromosome = re.sub("_P", "_T", chromosome)
                chromosome = re.sub("jgi\|", "", chromosome)
                #chromosome = re.sub("jgi\.p\|", "", chromosome)
                chromosome = re.sub("\|rna.*", "", chromosome)
                chromosome = re.sub("\|g\d*.*", "", chromosome)
                chromosome = re.sub("\|Chrsp.*", "", chromosome)
                chromosome = re.sub("\|Gphleg.*", "", chromosome)
                chromosome = re.sub("\|Mesvi.*", "", chromosome)
                chromosome = re.sub("\|ME.*", "", chromosome)
                chromosome = re.sub("\|CE.*", "", chromosome)
                chromosome = re.sub("\|estExt.*", "", chromosome)
                chromosome = re.sub("\|e_gw.*", "", chromosome)
                chromosome = re.sub("\|fgenesh.*", "", chromosome)
                chromosome = re.sub("\|MIX.*", "", chromosome)
                chromosome = re.sub("_\d_ORF\d", "", chromosome)
                #print >> sys.stderr, chromosome
                #if len(chromosome) > 20:
                #    print >> sys.stderr, chromosome
                if "Sh_" == chromosome[:3]:
                    chromosome = re.sub("_p", "_t", chromosome)
                if "GSMUA_A" == chromosome[:7]:
                    chromosome = re.sub("P", "T", chromosome)
                if "_cds_" in chromosome: #Pdactylifera, Sindicum, Thassleriana
                    chromosome = chromosome.split("_cds_")[1]
                    chromosome = re.sub("_\d$", "", chromosome)
                genomeDict.update({chromosome:[]})
            else:
                genomeDict[chromosome].append(line[0:-1])
        for chromosome in genomeDict:
            sequence = ''.join(genomeDict[chromosome])
            genomeDict[chromosome] = sequence
    return genomeDict

def loadPrimaryProteins(gffFile):
    noCDS = 0
    out = []
    skipFeatures = ["biological_region", "cDNA_match", "chromosome", "contig", "intron", "lnc_RNA", "miRNA", "ncRNA", "ncRNA_gene", "pre_miRNA", "region", "RNase_MRP_RNA", "rRNA", "scaffold", "snoRNA", "snRNA", "SRP_RNA", "start_codon", "stop_codon", "supercontig", "sequence_feature", "tmRNA", "tRNA", "transcription_end_site", "transcription_start_site", "transposable_element", "transposable_element_gene", "exon", "five_prime_UTR", "five_prime_utr", "three_prime_UTR", "three_prime_utr", "expressed_sequence_match", "UTR_3", "UTR_5", "intergenic", "inverted_repeat", "pseudogenic_tRNA", "repeat_region", "match", "match_part", "protein_match", "similarity", "polypeptide"]
    geneFeatures = ["gene", "pseudogene"]
    transcriptFeatures = ["mRNA", "pseudogenic_transcript", "transcript"]
    proteinFeatures = ["CDS"]
    isFirst = True
    fieldPatternExtended = re.compile("ID=.+?;|Name=.+?;|;Protein_Accession=.+?;$|protein_id=.+?;") # risky as this also extracts geneName=
    with myopen(gffFile, 'rb') as infile:
        isFirstWarning = True
        thisHasNoGene = False
        for line in infile:
            if line[:6] == "##FASTA":
                print >> sys.stderr, "found fasta entries, exiting loop"
                break
            if line[0] == "#":
                continue
            if len(line) < 30:
                continue
            line = line.strip()
            if line[-1] != ";":
                line += ";"
            try:
                (chrom, source, feature, start, end, score, strand, phase, rest) = line.split('\t')
            except ValueError:
                fields = line.split('\t')
                if len(fields) >= 9:
                    (chrom, source, feature, start, end, score, strand, phase, rest) = fields[:9]
                else:
                    print >> sys.stderr, fields
                    sys.exit(1)
            if source == "repeatmasker":
                continue
            if feature in skipFeatures:
                continue
            if "=" not in rest and ";" not in rest:
                rest = "ID="+rest+";"
            allMatches = fieldPatternExtended.findall(rest)
            rest = ';'.join(allMatches)
            # gene, store last
            if feature in geneFeatures:
                curGene = rest
                if isFirst:
                    isFirst = False
                else:
                    transcriptToSize.sort(reverse=True)
                    try:
                        out.append(transcriptToSize[0][1])
                    except IndexError:
                        # that's if there's no CDS
                        noCDS += 1
                transcriptToSize = []
                continue
            # transcript, just keep ID
            if feature in transcriptFeatures:
                curTranscript = rest
                isFirstCDS = True
                if thisHasNoGene:
                    try:
                        out.append(transcriptToSize[0][1])
                    except IndexError:
                        # that's if there's no CDS
                        noCDS += 1
                    transcriptToSize = [[0, ';'.join([curTranscript])]]
                elif "transcriptToSize" not in locals() and "transcriptToSize" not in globals():
                    print >> sys.stderr, "no gene ID!"
                    thisHasNoGene = True
                    transcriptToSize = [[0, ';'.join([curTranscript])]]
                else:
                    transcriptToSize.append([0, ';'.join([curGene, curTranscript])])
                continue
            if feature in proteinFeatures:
                #rest = re.sub("Parent=.*$", "", rest)
                try:
                    transcriptToSize[-1][0] += abs(int(end)-int(start))
                    if isFirstCDS:
                        isFirstCDS = False
                        if rest:
                            transcriptToSize[-1][1] = ';'.join([transcriptToSize[-1][1], rest])
                except IndexError:
                    if isFirstWarning:
                        isFirstWarning = False
                        print >> sys.stderr, "wrong sorting!"
                        print >> sys.stderr, line.strip()
        transcriptToSize.sort(reverse=True)
        try:
            out.append(transcriptToSize[0][1])
        except IndexError:
            # that's if there's no CDS
            noCDS += 1
    return out

allProteins = loadFasta(proteinFasta)
primaryProteins = loadPrimaryProteins(gffInfileName)
# clean up the rest matches
primaryProteinsClean = []
for curElem in primaryProteins:
    if ';' not in curElem:
        primaryProteinsClean.append(curElem)
        continue
    fields = [x.split('=') for x in curElem.split(';') if x]
    allMatches = ['='.join([key, val]) for key, val in fields if key in ["ID", "Name", "Protein_Accession", "protein_id"]]
    if len(allMatches) > 0:
        newEntry = ';'.join(allMatches)
        primaryProteinsClean.append(newEntry)
        #print >> sys.stderr, newEntry
    else:
        print >> sys.stderr, curElem

#fieldPattern = re.compile("ID=.+?;|;Protein_Accession=.+?$")
#fieldPatternExtended = re.compile("ID=.+?;|Name=.+?|;Protein_Accession=.+?$") # risky as this also extracts geneName=
#primaryProteinsClean = []
#for curElem in primaryProteins:
#    if ';' not in curElem:
#        primaryProteinsClean.append(curElem)
#        continue
#    #fields = curElem.split(';')
#    allMatches = fieldPattern.findall(curElem)
#    if len(allMatches) > 0:
#        newEntry = ''.join(allMatches)
#        primaryProteinsClean.append(newEntry)
#        #print >> sys.stderr, newEntry
#    else:
#        allMatches = fieldPatternExtended.findall(curElem)
#        if len(allMatches) > 0:
#            newEntry = ''.join(allMatches)
#            primaryProteinsClean.append(newEntry)
#        else:
#            print >> sys.stderr, curElem

#print >> sys.stderr, primaryProteinsClean[0]
primaryProteinsString = '|||'.join(primaryProteinsClean)
# plain, matches .1 to .10 
toKeep = [protName for protName in allProteins if protName in primaryProteinsString]
# with regex match
if len(toKeep) > len(primaryProteinsClean):
    print >> sys.stderr, "doing regex"
    toKeepClean = [protName for protName in toKeep if re.search("%s[^0-9a-zA-Z]" % protName, primaryProteinsString)]
else:
    toKeepClean = toKeep
#for protName in allProteins:
#    if protName not in primaryProteinsString:
#        print >> sys.stderr, protName
#for curName in allProteins:
#    print >> sys.stderr, curName
    
#alreadySeen = set([])
if len(toKeepClean) - len(primaryProteinsClean) < (-0.1*len(primaryProteinsClean)) or "Cmicranthum" in outfileName or "Dexilis" in outfileName or "Mbalbisiana" in outfileName or len(primaryProteinsClean) == 0:
    with myopen(outfileName, "w") as outfile:
        for curProt in allProteins:
            #curGene = curProt.split('.')[0]
            #if curGene in alreadySeen:
            #    print >> sys.stderr, "double: " + curProt
            print >> outfile, ''.join([">", curProt])
            print >> outfile, allProteins[curProt]
            #alreadySeen.add(curGene)
    print >> sys.stderr, "%s: wrote %d sequences as primary, found %d primary sequences in GFF, difference too large: %d (total in file: %d)" % (outfileName, len(allProteins), len(primaryProteinsClean), len(toKeepClean)-len(primaryProteinsClean), len(allProteins))
    if len(allProteins) == len(toKeepClean):
        print >> sys.stderr, "ok as all original are still in"
else:
    with myopen(outfileName, "w") as outfile:
        for curProt in toKeepClean:
            #curGene = curProt.split('.')[0]
            #if curGene in alreadySeen:
            #    print >> sys.stderr, "double: " + curProt
            print >> outfile, ''.join([">", curProt])
            print >> outfile, allProteins[curProt]
            #alreadySeen.add(curGene)
    print >> sys.stderr, "%s: wrote %d sequences but found %d primary sequences in GFF, difference: %d (total in file: %d)" % (outfileName, len(toKeepClean), len(primaryProteinsClean), len(toKeepClean)-len(primaryProteinsClean), len(allProteins))
    if len(allProteins) == len(toKeepClean):
        print >> sys.stderr, "ok as all original are still in"












