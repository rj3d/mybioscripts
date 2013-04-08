import os
import time
import glob
import tempfile as temp
import threading
import subprocess as sp
import sys

def parseFasta(fasta):
    f = open(fasta)
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
                sequences.append( (cur_title, ''.join(cur_sequence) ) )
                del cur_sequence[:]
            cur_title = cur_line[1:].strip()
        else: #If still in current sequence
            cur_sequence.append(cur_line.strip())
    if not cur_sequence == '': #Append infor for final sequence
        sequences.append( (cur_title, ''.join(cur_sequence) ) )
    return sequences

def getCmds(fasta, seqs, outdir):
    stack = []
    i = 0
    for seq_name, seq in seqs:
        i += 1
        pat_id = seq_name.split()[0]
        stack.append(['patmatmotifs', '-sequence', ':'.join([fasta,pat_id]), '-outfile', os.path.join(outdir, str(i)), '-full'])
    return stack

#Worker to run a command and wait for it to terminate.
def worker(cmd):
    proc = sp.Popen(cmd, bufsize = -1, cwd = tempdir)
    proc.wait()
    return

def runCmds(max_processes, cmds):
    max_threads = max_processes  + 1 #Need to account for this thread, the main thread
    #Run through all commands and fill processes dict with running Popen processes
    while len(cmds):
        cmd = cmds.pop()
        t = threading.Thread(target=worker, args=(cmd,))
        t.start()
        #Check if we need to start new threads every 60 seconds.
        while threading.active_count() == max_threads:
            time.sleep(1)

def parsePMM(pmm):
    f = open(pmm)
    #Get to the start of the summary section
    name = ''
    length = -1
    count = 0
    while count < 2:
        cur_line = f.readline()
        #Get the name of the sequence
        if 'Sequence:' in cur_line:
            splat = cur_line.split()
            name = splat[2]
            length = int(splat[6])
        if '#=======================================' in cur_line:
            count += 1
    #Now in the summary section, you have 'Length = ' followed by 'Motif = ' a few lines later.
    motifs = []
    cur_length = -1
    while True:
        cur_line = f.readline()
        if '#---------------------------------------' in cur_line:
            break
        if 'Length =' in cur_line:
            splat = cur_line.split()
            cur_length = int(splat[2])
        if 'Motif =' in cur_line:
            splat = cur_line.split()
            motifs.append((splat[2], cur_length))
    return name, length, motifs
        
def parseOutputs(tempdir):
    info = []
    listing = os.listdir(tempdir)
    for infile in listing:
        name, length, motifs = parsePMM(os.path.join(tempdir, infile))
        info.append((name, length, infile, motifs))
    return info

def writeOutputs(info, outdir):
    summary = open(os.path.join(outdir, 'summary.tsv'),'w')
    summary.write('Name\tLength\tScore\tRaw_File\tMotifs\n')
    stats = {}
    lens = {}
    all_len = 0
    num_seqs = 0
    #Iterate through the sequences and write the summary file
    for name, seq_len, raw, motifs in info:
        all_len += seq_len
        num_seqs += 1
        summary.write('\t'.join([name, str(seq_len), '']))
        score = 0.0
        #Iterate through the motifs for the current sequence
        for motif, motif_len in motifs:
            score += motif_len
            #Get motif distribution info
            if not motif in stats:
                stats[motif] = 0
                lens[motif] = motif_len
            stats[motif] += 1
        summary.write('\t'.join([str(score/seq_len), raw, '']))
        summary.write(', '.join([motif for motif, motif_len in motifs]))
        summary.write('\n')
    #Write distribution info file
    dist = open(os.path.join(outdir, 'distribution.tsv'),'w')
    dist.write('Stats: \t')
    dist.write('Total bp: ' + str(all_len) + '\t')
    dist.write('Num seqs: ' + str(num_seqs) + '\t')
    dist.write('\n')
    dist.write('Motif\tCount\tLength\n')
    for motif in stats.iterkeys():
        dist.write('\t'.join([motif, str(stats[motif]), str(lens[motif]) ]))
        dist.write('\n')
        
if __name__ == "__main__":
    if len(sys.argv) == 3:
        fasta = os.path.abspath(sys.argv[1])
        outdir = os.path.abspath(sys.argv[2])
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        tempdir = os.path.join(outdir, 'raw')
        if not os.path.exists(tempdir):
            os.mkdir(tempdir)
        seqs = parseFasta(fasta)
        cmds = getCmds(fasta, seqs, tempdir)
        runCmds(15, cmds)
        #Wait for the last threads to finish before exiting.
        while threading.active_count() > 1 :
            time.sleep(1)
        info = parseOutputs(tempdir)
        writeOutputs(info, outdir)
    else:
        print "\nUsage:"
        print "python parallelPMM.py input.fasta ./outdir/"
        print "Runs patmatmotifs on every sequence in the fasta file."
        print "Generates outdir if it doesn't exist, and then generates summary and distribution files.\n"

