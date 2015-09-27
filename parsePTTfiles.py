#! usr/bin/python2.7

"""
 Author: Jagannath
 Date: 25th Sept 2015
 Last modified : 25th Sept 2015
 
 The file parses the .ptt file to extract information about the CDS (protein coding genes). Locus_tag is used as the identifier although PID (protein ID) is also provided. The dictionaries created are - 
 (1) { locus_tag : [accNbr, rank, geneName, PID (protein ID), (startPos,endPos), COG (COG number with Code)] }
 (2) { COG (w/ code) : ( [locus_tag1, locus_tag2, ...],Count of locus_tags ) }
 (3) { COG(A)-COG(B) : count } # Count of number of times when COG(A) before COG(B); A --> B on the same strand
 [! In future, maybe can think of A --> B (but A and B are on different strands)
"""

import os
from collections import defaultdict
import fnmatch
import pickle

class PTTfile:

    def __init__(self,pttFile):
        self.pttFile = pttFile
        self.accNbr = os.path.split(pttFile)[-1].split('.')[0]
    
    def openPttFile(self):
        f = self.pttFile
        ifile = open(f,'r')
        lines = ifile.readlines()[3:] #Skipping the header
        ifile.close()
        return lines
        
    def getPos(self,loc,strand):
        if strand is '+':
            startPos, endPos = map(int,loc.split('..'))
        else:
            endPos, startPos = map(int,loc.split('..'))
        return (startPos,endPos)

    def makeLinePairs(self,lines):
        # Makes pairs of adjacent genes; By default all genomes are considered circular. So the last pair in the lines 
        # is the last gene --> first gene
        # Returns list of [(line0,line1,rank0),(line1,line2,rank1),...(line4175,line0,rank4175)]
        nbrGenes = len(lines)
        line_first = lines[0]
        line_last = lines[-1]
        lines_shift = lines[1:]
        lines_shift.append(line_first)
        linePairs = zip(lines,lines_shift,range(nbrGenes))
        return linePairs

    def getCDSinfo(self,line,rank):
        # Column information
        # [0]Location	[1]Strand	[2]Length	[3]PID	[4]Gene	[5]Synonym	[6]Code	[7]COG	[8]Product
        (loc, strand, pid, gene, locus_tag, COG) = (line.split('\t')[i] for i in [0,1,3,4,5,7])
	if ',' in COG: # There are some wierd .ptt files wherein the COG is written as COG1,COG1
	    COG = COG.split(',')[0]
        #pos = getPos(loc,strand)
        pos = map(int,loc.split('..'))
        cds_info = [self.accNbr, rank, locus_tag, gene, pid, pos, strand, COG]
        return cds_info

    def getCOGpair(self,geneA,geneB):
        # Returns a tuple of COGpairs; If orientation is '-' and '-', then the cogpairs are reversed accordingly
        cogPair = ('-','-')
        strA = geneA[-2]
        strB = geneB[-2]
        if strA == strB == '+':
            cogPair = (geneA[-1],geneB[-1])
        elif strA == strB == '-':
            cogPair = (geneB[-1],geneA[-1])
        else:
            cogPair = ('-','-')
        return cogPair

    def populateDict(self):
        lines = self.openPttFile()
        linePairs = self.makeLinePairs(lines)
        for pair in linePairs:
            cogPair = ('-','-')
            lineA,lineB,rankA = pair
            geneA = self.getCDSinfo(lineA,rankA)
            geneB = self.getCDSinfo(lineB,rankA+1)
            cogPair = self.getCOGpair(geneA,geneB)
            cogPair_count_dict[cogPair]+=1
            cog_locusTagList_dict[geneA[-1]].append(geneA[2])
    
# Defining dictionary
cogPair_count_dict = defaultdict(int)
cog_locusTagList_dict = defaultdict(list)

sourceDir = '/home/jaggu/research/allGenomePttFiles'
pklPath = '/home/jaggu/research/projectFiles/operons/pklFiles'

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and
    below supplied root directory.'''
    allFiles = []
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            allFiles.append(os.path.join(path,filename))
    return allFiles

def savePkl(db,pklFname):
    f = os.path.join(pklPath,pklFname)
    pickle.dump(db,open(f,'w'))
    return 

def main(sourceDir):
    allPttFiles = locate('*.ptt',sourceDir)
    for pttFile in allPttFiles:
        print "Parsing :  %s file ..."%(pttFile)
        genome = PTTfile(pttFile)
        genome.populateDict()
    savePkl(cogPair_count_dict,'cogPair_count.dict.pkl')
    savePkl(cog_locusTagList_dict,'cog_locusTag.dict.pkl')
        
def test_case():
    orgName = 'Bacillus_subtilis_168_uid57675'
    pttFile = 'NC_000964.ptt'
    test = PTTfile(os.path.join(sourceDir,orgName,pttFile))
    test.populateDict()
    print cogPair_count_dict

main(sourceDir)

