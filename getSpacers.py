import sys
import re
import itertools
import string

from collections import defaultdict

from fastarator import fasta

#Default values assuming Pyrogenes Cas9 PAM
#No mutations are tolerated in the first 11bp. Cong et. al, Science, January 2013
#~7 bp of mutations are required in the 5' end to ablate cutting if the core sequences are identical. Jinek et. al, Science, August 2012. 
class spacerParser:
    def __init__(self, pam_size = 3, pam = "GG", core_size = 11, end_size = 9, min_mismatch = 7):
        self.pam_size = pam_size
        self.pam = pam
        self.pam_rc = self.__revComp(pam)
        self.core_size = core_size
        self.end_size = end_size
        self.min_mismatch = min_mismatch

    def __revComp(self, sequence):
        complement = string.maketrans('ATGCN', 'TACGN')
        return sequence.translate(complement)[::-1]

    def parseSpacers(self, seq, title):
        #Offsets from reported regex start site for the left/right bounds of the targeting sequence for forward and reverse strand.
        f_left = self.core_size + self.end_size + 1
        f_right = 1
        r_left = self.pam_size + self.core_size + self.end_size - 1
        r_right = 0
        #Generates a list of tuples structured as follows (spacer_start, spacer_end, spacer_sequence) for both the forward and reverse strands.
        spacers = [s for s in itertools.chain(
                  [ (m.start()-f_left, m.start()+f_right, title, seq[m.start()-f_left:m.start()+f_right+1]) for m in re.finditer('(?='+ self.pam +')', seq)],
                  [ (m.start()+r_left, m.start()+r_right, title, self.__revComp(seq[m.start():m.start()+r_left+1])) for m in re.finditer('(?='+ self.pam_rc + ')', seq)])
                  ]
        return spacers

#For defaultdict in compareEnds
def defaultFact():
    return -1

def compareEnds(spacers, min_mismatch):
    filt = []
    min_scores = defaultdict(defaultFact)
    #Compare all unique combinations of spacers and save the lowest score (closest spacer)
    for i, j in itertools.combinations(spacers, 2):
        end_i = i['end']
        end_j = j['end']
        score = sum( [ 0 if x==y else 1 for x,y in itertools.izip(end_i,end_j)] ) #Count the number of mismatches
        if score < min_scores[end_i]: min_scores[end_i] = score
        if score < min_scores[end_j]: min_scores[end_j] = score
    #Filter the scored spacers based on min_mismatch
    for spacer in spacers:
        min_score = min_scores[spacer['end']]
        if min_score >= min_mismatch:
            spacer['min_score'] = min_score
            filt.append(spacer)
    return filt

#Given a dict of spacers keyed on core sequences, compares and filters all spacers with the same core.
def filterSpacerDict(spacer_dict, min_mismatch):
    final_spacers = []
    for core in spacer_dict.iterkeys():
        cur_spacers = spacer_dict[core]
        if len(cur_spacers) == 1: #If only one spacer for this core, no off-targeting
            cur_spacer = cur_spacers[0]
            final_spacers.append(cur_spacer)
        else: #Multiple spacers for the core, have to compare and look for off-targeting
            filt_spacers = compareEnds(cur_spacers, min_mismatch)
            final_spacers.extend(filt_spacers)
    return final_spacers

#Given a set of spacers and the parser used to generate them, produces a dict of spacers keyed on the core targeting sequence. 
def getSpacerDict(spacers, parser):
    spacer_dict = defaultdict(list)
    #Constants for the spacer based on the info in the parser
    spacer_core_start = parser.pam_size + parser.core_size
    spacer_len = spacer_core_start + parser.end_size
    for start, stop, title, spacer in spacers:
        #end is left most portion of spacer. PAM is right most portion. Core is between PAM and end.
        end = spacer[0:spacer_len-spacer_core_start]
        core = spacer[spacer_len-spacer_core_start:spacer_len-parser.pam_size]
        spacer_dict[core].append(
        {'end':end, 'start':start, 'stop':stop, 'title':title, 'spacer':spacer, 'min_score':-1}
        )    
    return spacer_dict    

#Finds all non-off targeting spacers in a set of sequences
def findSpacers(sequences, parser):
    spacers = []
    print "Generating spacers..."
    for title, sequence in sequences:
        spacers.extend(parser.parseSpacers(sequence, title))
    print "Generated " + str(len(spacers)) + " spacers."
    spacer_dict = getSpacerDict(spacers, parser)
    del spacers
    print "Collapsed into " + str(len(spacer_dict)) + " core sequences."
    final_spacers = filterSpacerDict(spacer_dict, parser.min_mismatch)
    del spacer_dict
    print "Found + " + str(len(final_spacers)) + " spacers without predicted off-targets."
    return final_spacers

def writeSpacers(spacers, outfile):
    f = open(outfile, 'w')
    f.write('Spacer\tSource\tStart\tStop\tMin_score\n')
    for spacer in spacers:
        f.write(spacer['spacer'])
        f.write('\t')
        f.write(spacer['title'])
        f.write('\t')
        f.write(str(spacer['start']))
        f.write('\t')
        f.write(str(spacer['stop']))
        f.write('\t')
        f.write(str(spacer['Min_score']))
        f.write('\n')
    f.close()
        
if __name__ == "__main__":
    if len(sys.argv) == 3:
        parser = spacerParser()#Can change the default options if desired.
        print "Parsing sequences... "
        sequences = fasta(sys.argv[1])
        spacers = findSpacers(sequences, parser)
        writeSpacers(spacers, sys.argv[2])
    else:
        print "Usage: python getSpacers.py input.fasta output.tsv"


