from collections import Iterator

# Iterator object for fasta file
# Returns tuple containing title, sequence.
class fasta(Iterator):

    def __init__(self, infile):
        self.seqs = self.__parseFasta(infile)
        self.num_seqs = len(self.seqs)
        self.current = -1

    def contents(self):
        return self.seqs

    def next(self):
        if self.current < self.num_seqs-1:
            self.current += 1
            return self.seqs[self.current]
        else:
            raise StopIteration

    #Parses a fasta file with either one line or multiple lines per sequence.
    def __parseFasta(self, infile):
        f = open(infile)
        sequences = []
        titles = []
        cur_sequence = []
        cur_title = ''
        while True:#Read through file
            cur_line = f.readline()
            if cur_line == '': #If EOF
                break
            elif '>' in cur_line: #If next sequence reached
                if not len(cur_sequence) == 0:
                    sequences.append( (cur_title, ''.join(cur_sequence).upper() ) )
                    del cur_sequence[:]
                cur_title = cur_line[1:].strip()
            else: #If still in current sequence
                cur_sequence.append(cur_line.strip().replace('-',''))
        if not cur_sequence == '': #Append info for final sequence
            sequences.append( (cur_title, ''.join(cur_sequence).upper() ) )
        return sequences

# Use example
if __name__ == "__main__":
    import sys
    ff = fasta(sys.argv[1])
    print "Successfully read " + str(ff.num_seqs) + " sequences."
    for title, seq in iter(ff):
        print "Iterating over sequence " + title
                
