import re
import string
import sys
from fastarator import fasta
from cachedwriter import CachedWriter

#Finds all of the substrings in the passed in substrings list in the passed in string.
def findSubstrings(string, substrings):
    indices = []
    for substring in substrings:
        indices.extend([m.start() for m in re.finditer(substring, string)])
    return sorted(indices, reverse=True)

#Takes the first stop codon then looks at all starts prior to that stop.
#If any are in frame, creates an ORF with that start and stop pair.
def findORFs(sequence, min_len, max_len):
    starts = findSubstrings(sequence, ['CTG', 'ATG'])
    stops =  findSubstrings(sequence, ['TAG', 'TAA', 'TGA'])
    indices= []
    if len(starts) > 0 and len(stops) > 0:
        for cur_stop in stops[::-1]: #For a given stop codon, assume that all start codons in frame are potential smORFs.
            while  len(starts) > 0 and starts[-1] > cur_stop: #Remove all start codons that are no longer needed
                starts.pop()
            for cur_start in starts[::-1]:
                if cur_start > cur_stop: #If next start is past the current stop.
                    break
                if cur_stop%3 == cur_start%3:#If in frame, remove that start.
                    if min_len < (cur_stop - cur_start) < max_len: #If correct length, add to list.
                        indices.append((cur_start, cur_stop))
                    starts.remove(cur_start)            
    return indices

def getSubstrings(sequence, indices, up, down):
    substrings = []
    for start, stop in indices:
        cur_seq = ''.join([sequence[start-up:start].lower(), sequence[start:stop+3].upper(), sequence[stop+3:stop+down].lower()])
        substrings.append( cur_seq )
    return substrings
    
def reverseComplement(sequence):
    complement = string.maketrans('ATGCN', 'TACGN')
    return sequence.translate(complement)[::-1]

def getORFs(sequence, min_len, max_len, upstream, downstream):
    fi = findORFs(sequence, min_len, max_len)
    orfs = getSubstrings(sequence, fi, upstream, downstream)
    rc = reverseComplement(sequence)
    ri = findORFs(rc, min_len, max_len)
    orfs.extend(getSubstrings(rc, ri, upstream, downstream))
    return orfs

#NOTE: Want small cutoff ~8-12aa.
#Max length is 210 nuc, or 70 aa.
#244k - 488k oligos max for synthesis
if __name__ == "__main__":
    #Parse fasta. Iterate over sequences, then getORFs.
    sequences = fasta(sys.argv[1])
    #Magic numbers
    upstream = 0
    downstream = 0
    min_aa = 8
    min_nuc = min_aa * 3
    max_aa = 70
    max_nuc = max_aa * 3
    #Init variables
    cw = CachedWriter(sys.argv[2])
    orf_num = 0
    orfs = {}
    #GO!
    for title, sequence in sequences:
        cur_num = 0
        for orf in getORFs(sequence, min_nuc, max_nuc, upstream, downstream):
            if not orf in orfs:#filter out duplicates
                cur_num += 1
                orf_num += 1
                orfs[orf] = orf_num
                cw.write(''.join(['>smORF_', str(orf_num), ' ORF_', str(cur_num), ' Parent_', title]))
                cw.write('\n')
                cw.write(orf)
                cw.write('\n')
    cw.flush()

    
