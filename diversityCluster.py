# -*- coding: utf-8 -*-
"""
Created on Thu May 19 11:39:58 2016

@author: Ben Demaree
"""

### Misc. imports
from __future__ import division
import numpy as np
import os
import sqlite3
from random import randint, shuffle
import math
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from multiprocessing import Process, Queue
import argparse
import os.path

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

### Pandas imports
import pandas as pd
from sqlalchemy import create_engine

####################################################################################

def encodeSeq (seq):
# Encodes a DNA sequence string to a bitstring representation (list format)

    # The 4-bit encoding scheme:
    encoding = {'A': '0001', 'C': '0010', 'G': '0100', 'T': '1000'}
    bases = ['A', 'C', 'G', 'T']
    
    # Map bases to encoding dictionary:
    
    encodedSeq = ''     # Encoded sequence

    for base in seq:
        
        # Map ambiguous bases (N) to random nucleotide
        if base == 'N':
            encodedSeq += encoding[bases[randint(0, 3)]]
        
        # Map known bases to corresponding bitstring:
        else:
            encodedSeq += encoding[base]      
    
    return encodedSeq
    
def extractBargroup (barcode, table):
# For a given barcode group, extracts all sequences (including both mate-pairs)
# from the sqlite DB, encodes them, and stores them as a list of lists.

    encodedSeqs = []        # List for storing encoded bargroup sequences
    originalSeqs = []       # List for storing original sequences

    # For reverse-complementing sequences:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'n': 'n', 'N': 'N'}
    
    if not SE:
        seqs = pd.read_sql_query('select seq1,seq2 from "%s" where barcode="%s"' % (table, barcode), db)
        seq1 = seqs['seq1']
        
        # Reverse-complement option:
        if RC:
            seq2 = "".join(complement.get(base, base) for base in reversed(seqs['seq2']))        
        else:
            seq2 = seqs['seq2']
    
    else:
        seqs = pd.read_sql_query('select seq1 from "%s" where barcode="%s"' % (table, barcode), db)
        seq1 = seqs['seq1']
    
    for i in range(len(seq1)):
        
        # Encode sequence strings and add to list:
        encodedSeqs += [encodeSeq(seq1[i])]
        if not SE:
            encodedSeqs += [encodeSeq(seq2[i])]
        
        # Also store original sequences in same order:
        originalSeqs += [seq1[i]]
        if not SE:
            originalSeqs += [seq2[i]]

    return encodedSeqs, originalSeqs
    
def removeN (seq):
# If a sequence contains an ambiguous nucleotide (N), replaces it with a random 
# nucleotide (A, T, C, G)

    bases = ['A', 'C', 'G', 'T']

    if 'N' not in seq:
        return seq
        
    else:
        # Find indices of Ns:
        Ns = [i for i in range(len(seq)) if seq[i] == 'N']
        
        # Replace Ns with random base:
        for j in Ns: 
            seq = seq[:j] + bases[randint(0, 3)] + seq[j+1:]
            
        return seq
        
def appendZero (bitstring):
# Appends 0s to the end of an encoded sequence up to a specified target length

    lenDiff = targetLen * 4 - len(bitstring)
    
    bitstring += '0' * lenDiff
    
    return bitstring
    
def linkageMat (encodedSeqs):
# Receives encoded sequences, appends 0s, computes pairwise distances, and
# generates a linkage matrix.

    X = []      # Array for holding sequences

    for seq in encodedSeqs:
        seqZ = appendZero(seq)
        X += list(seqZ)
    
    # Reshape input vectors to correct size:
    X = np.array(X, dtype = int)
    X = np.reshape(X, (len(encodedSeqs), targetLen * 4))
    
    # Generate linkage matrix:
    Z = linkage(X, method = 'average', metric = 'hamming')
    
    return Z
    
def plotDendrogram (Z, barcode, encodedSeqs):
# Plots a dendrogram from the linkage matrix for the specified bargroup.

    plt.figure()
    plt.title('Hierarchical Clustering Dendrogram (Truncated)\nBarcode = %s - # Reads = %d' 
    % (barcode, len(encodedSeqs)))
    
    plt.xlabel('Sample Index or (Cluster Size)')
    plt.ylabel('Hamming Distance')
    dendrogram(
        Z,
        truncate_mode = 'lastp',    # Show only the last p merged clusters
        p = 12,                     # Show only the last p merged clusters
        leaf_rotation = 90.,
        leaf_font_size = 12.,
        show_contracted = True,     # To get a distribution impression in truncated branches
    )
    plt.gcf().subplots_adjust(bottom = 0.15)
    plt.savefig('dendrogram_%s.pdf' % barcode)
    
def outputClusters (Z, barcode, originalSeqs, textOutput = True):
# Outputs a text file containing the clustered sequences

    # Flatten the clusters:
    T = fcluster(Z, hamming, 'distance')
    # Judging by the cluster lists, d = 0.05 is a conservative estimate for determining 
    # read diversity.

    numClusters = max(T)
    
    if textOutput:
        
        outFile = open('clusters_%s.txt' % barcode, 'w')
        outFile.write('ID\tSequence\n')
        
        for i in range(1, numClusters + 1):
            
            # Find sequences in ith cluster:
            ind = [j for j in range(len(T)) if T[j] == i]
            
            outFile.write('\n----- Cluster %d (%d reads) -----\n' % (i, len(ind)))
    
            for k in ind:
                outFile.write('%d\t%s\n' % (k, originalSeqs[k]))
            
        outFile.close()
    
    return numClusters
    
def randomRead (length):
# Generates a random read of specified length

    bases = ['A', 'C', 'G', 'T']
    
    read = ''    
    
    for i in range(length):
        
        read += bases[randint(0, 3)]
        
    return read
    
def cluster (barcode):
# For a given barcode group, clusters reads and outputs read statistics

    originalSeqs = []
    encodedSeqs = []
    
    encodedSeqs, originalSeqs = extractBargroup(barcode, table)
    
    origLen = len(originalSeqs)         # Original number of reads
    
    # Remove exact duplicates:
    originalSeqs = list(set(originalSeqs))  
    encodedSeqs = list(set(encodedSeqs))
    
    uniqLen = len(originalSeqs)         # Number of unique reads
    
    Z = linkageMat(encodedSeqs)         # Compute linkage matrix
    numClusters = outputClusters(Z, barcode, originalSeqs, False)   # Number of clusters
    
    output = '%d_%d_%d' % (origLen, uniqLen, numClusters)
    
    return output
    
def mpCluster (barcodes, nprocs):
    
    def worker (barcodes, out_q):
        """ The worker function, invoked in a process. 'barcodes' is a
            list of barcodes to cluster. The results are placed in
            a dictionary that's pushed to a queue.
        """
        
        outDict = {}
        
        for bar in barcodes:
            outDict[bar] = cluster(bar)
            
        out_q.put(outDict)

    # Each process will get random subset of barcodes:
    out_q = Queue()
    chunkSize = int(math.ceil(len(barcodes) / float(nprocs)))
    procs = []
    
    shuffle(barcodes)       # Shuffle barcodes
    
    if os.path.isfile('libraryDiversity_%s.txt' % expName):
        print 'Output text file exists! Please rename or delete the existing file.'
        raise SystemExit
    
    else:   
        # Prepare to write clustering data to a text file:
        outFile = open('libraryDiversity_%s.txt' % expName, 'w')
        outFile.write('Barcode\tTotal Reads\tUnique Reads\tRead Clusters\n')

    for i in range(nprocs):
        p = Process(
                target=worker,
                args=(barcodes[chunkSize * i:chunkSize * (i + 1)],
                      out_q))
        procs.append(p)
        p.start()

    # Collect all results into a single result dict:
    resultDict = {}
    
    for i in range(nprocs):
        resultDict.update(out_q.get())

    # Wait for all worker processes to finish:
    for p in procs:
        p.join()
    
    for key in resultDict:
        clusterData = resultDict[key].split('_')
        origLen = int(clusterData[0])
        uniqLen = int(clusterData[1])
        numClusters = int(clusterData[2])
        outFile.write('%s\t%d\t%d\t%d\n' % (key, origLen, uniqLen, numClusters))

    outFile.close()    
    
    return resultDict
    
def insertDiversity (barcodes, diversity, table):
# Inserts read diversity information into the sqlite database.

    sql = sqlite3.connect(sqliteDB, isolation_level = 'Exclusive')
    s = sql.cursor()    # Create cursor
    # These settings allow DB to be stored in memory and written periodically to the disk: 
    s.execute('PRAGMA synchronous = 0')
    s.execute('PRAGMA journal_mode = OFF')
    
    # Create index on barcodes (if not already done):
    try:
        s.execute('create index barIndex on "%s" (barcode)' % table)       
    except sqlite3.Error:
        pass

    # Add column for diversity (if not already done):
    try:
        s.execute('alter table "%s" add column diversity real' % table)        
    except sqlite3.Error:
        pass
            
    # Insert diversity information:
    for i in range(len(barcodes)):            
        s.execute('update "%s" set diversity="%s" where barcode="%s"' %
        (table, diversity[i], barcodes[i]))
        
        if (i % 1e3 == 0 and i != 0):
            sql.commit()        # Write reads to persistent DB
            print "%d diversity scores exported to database so far..." % i
                
    sql.commit()        # Commit db changes to disk one last time
    sql.close()         # Close the db
    
    print '%d barcode group diversities inserted into database!' % len(barcodes)
    

if __name__ == "__main__":
    
    ### Setup argument parser ###
     
    parser = argparse.ArgumentParser(description="""
    Clusters the reads within barcode groups to obtain a measure of sequencing
    diversity. Diversity scores are added to the specified SQLite database.
    
    Alternatively, the cluster file and dendrogram can be produced for a single 
    specified barcode.
    
    When the sequences in the SQLite DB are from different strands, select the 
    --RC flag to reverse-complement one of the mates before clustering.
    
    ***Dependencies***
    sqlite3
    
    @author: Ben Demaree
    """, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('db', metavar='<db_filename>', type=str,
                        help='SQLite input database filename')
                        
    parser.add_argument('read_length', type=int, metavar='length',
                        help='maximum read length (e.g. 50, 75, 150, etc...)')
                        
    parser.add_argument('-hd', type=float, metavar='hamming_distance', default=0.05,
                        help='Hamming distance for clustering reads [default: 0.05]')
                        
    parser.add_argument('-b', metavar='barcode', type=str, default='none',
                        help='barcode to use for cluster file and dendrogram plot')
                        
    parser.add_argument('-t', metavar='table', type=str, default='reads',
                        help='table in database to analyze [default: reads]')
                        
    parser.add_argument('-nprocs', type=int, metavar='N', default=4,
                        help='number of processes to run [default: 4]')             
                            
    parser.add_argument('--SE', action='store_true', default=False,
                        help='flag for single-end read analysis') 
                        
    parser.add_argument('--RC', action='store_true', default=False,
                        help='reverse-complement one of the paired-end reads before clustering') 
    
    args = parser.parse_args()      # Parse arguments
        
    sqliteDB = args.db              # DB filename
    targetLen = args.read_length    # Size of all sequences (shorter sequences padded)
    table = args.t                  # Table in DB to analyze
    nprocs = args.nprocs            # Number of processes to run
    bar = args.b                    # Barcode for outputting cluster file and dendrogram
    hamming = args.hd               # Hamming distance threshold for clustering reads
    SE = args.SE                    # Option to analyze a single-end read
    RC = args.RC                    # Option to RC one mate from a paired-end read
    
    expName = sqliteDB.split('.')[0]    # Universal filename label
    
    db = create_engine('sqlite:///' + sqliteDB)     # Create DB interface
    
####################################################################################  
# Option 1: Generate cluster file and dendrogram for a specified barcode
#################################################################################### 
    
    if bar != 'none':
        
        encodedSeqs, originalSeqs = extractBargroup(bar, table)
        print '%s: this barcode group contains %d reads.' % (bar, len(encodedSeqs))
        
        # Compute linkage matrix and plot dendrogram:
        Z = linkageMat(encodedSeqs)
        plotDendrogram(Z, bar, encodedSeqs)
        print 'Dendrogram plotted!'
        
        # Output cluster file:
        numClusters = outputClusters(Z, bar, originalSeqs)
        print 'Clusters exported to file!'
        
####################################################################################  
# Option 2: Cluster reads in all bargroups and export to text file and database
####################################################################################         
    
    else:
        
        print 'Initializing %d processes...' % nprocs
        
        print 'Beginning multiprocessor clustering...'
        
        bars = pd.read_sql_query('select barcode from "%s" group by barcode' % table, db)
        barcodes = bars['barcode']
        resultDict = mpCluster(barcodes, nprocs)
        
        print '%d barcode groups clustered! Text file saved.' % len(resultDict)
    
        # Insert diversity scores into database:
    
        print 'Inserting diversity scores into database...'
        
        libDiversity = open('libraryDiversity_%s.txt' % expName, 'r')
    
        barcodes = []
        diversity = []
        
        libDiversity.next()         # Skip header line
    
        for line in libDiversity:
            
            barcodes += [line.split('\t')[0]]
            diversity += [int(line.split('\t')[3]) / int(line.split('\t')[1])]
            
        libDiversity.close()
    
        insertDiversity(barcodes, diversity, table)

####################################################################################  
# Supplemental: Clustering random reads
####################################################################################    
      
#    ##### Clustering Randomly-Generated Reads #####
#        
#    numReads = 500
#    readLen = 200
#    
#    originalSeqs = []
#    encodedSeqs = []
#    
#    for i in range(numReads):
#        read = randomRead(readLen)
#        originalSeqs += [read]
#        encodedSeqs += [encodeSeq(read)]
#        
#    # Compute linkage matrix and plot dendrogram:
#    Z = linkageMat(encodedSeqs)
#    plotDendrogram(Z, 'random', encodedSeqs)
#    
#    # Output cluster file:
#    numClusters = outputClusters(Z, 'random', originalSeqs)
#    
#    # Result: 500 different "clusters"!!!





        
        
    
    
    
    

    