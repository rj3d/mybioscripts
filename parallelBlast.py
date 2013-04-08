import os
import time
import glob
import tempfile as temp
import threading
import subprocess as sp
import sys

#Worker to run a command and wait for it to terminate.
def worker(cmd):
    infile = cmd[6]
    tempdir =  temp.mkdtemp()
    print "STATUS: Worker " + str(threading.currentThread()) + " starting BLAST on file " + infile
    proc = sp.Popen(cmd, bufsize = -1, cwd = tempdir)
    proc.wait()
    print "STATUS: Worker " + str(threading.currentThread()) + " completed BLAST on file " + infile
    return
    
#Prepare the commands to run in parallel
inpath = sys.argv[2]
outpath = sys.argv[3]
database = sys.argv[1]
stack = []
for infile in glob.glob( os.path.join(inpath, '*.fasta') ):
    basename = os.path.basename(infile)
    filename = basename #filename = basename.split('.')[0]
    outfile = ''.join([outpath,filename])
    stack.append(['blastall', '-p', 'blastx', '-d', database, '-i', infile, '-e', '.001', '-o', outfile, '-m', '8', '-a', '2', '-W', '3'])
print "STATUS: Initialized stack with " + str(len(stack)) + " files."

#Start the commands in new threads.
max_processes = 64
max_threads = max_processes  + 1 #Need to account for this thread, the main thread
#Run through all commands and fill processes dict with running Popen processes
while len(stack):
    cmd = stack.pop()
    t = threading.Thread(target=worker, args=(cmd,))
    t.start()
    #Check if we need to start new threads every 60 seconds.
    while threading.active_count() == max_threads:
        time.sleep(60)

#Wait for the last threads to finish before exiting.
while threading.active_count() > 1 :
    time.sleep(60)
        
print "STATUS: Parallel blast complete"

