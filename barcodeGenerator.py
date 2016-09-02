# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 09:27:57 2016

@author: Ben Demaree
"""

'''
Generates simulated groups of barcoded reads from reference genomes in FASTA format.
'''

import random

def fasta2bargroups(fasta, readLen, barSize, numBargroups):
# Extracts a selected number of reads from randomly selected regions in a FASTA
# sequence file. Outputs the reads as a bargroup FASTA.

    f = open(fastaDir + fasta)  # Open FASTA
    genomes = f.read()          # Save FASTA in string
    numBases = len(genomes)     # Number of chars in FASTA (incl. headers and sequences)
    
    bases = 'ATCGatcg'          # Allowable bases
    
    outFile = 'bargroups_' + fasta      # Output bargroup FASTA
        
    genomeSeqs = []     # List of lists containing sequences for each synthetic bargroup
    
    for i in range(numBargroups):
            
        seqs = []       # Extracted sequences
        numSeqs = 0     # Number of extracted sequences
        
        while numSeqs < barSize:   # Until we get enough valid sequences
            
            goodSeq = True      # Flag for good seqs (not containing header chars)
            
            index = random.randint(0, numBases)         # Random index to extract
            seq = genomes[index:(index + readLen)]
            seq = seq.replace('\n', '')                 # Remove newline characters
            
            for b in seq:
                if b not in bases:
                    goodSeq = False     # Seq contains header, throw it out
                    break
                
            if goodSeq:
                seqs += [seq]
                numSeqs += 1
         
        genomeSeqs += [seqs]
        
    # Write barcode groups to FASTA file:
        
    out = open(outFile, 'w')
    
    header = fasta.split('.')[0]
    
    bargroupCount = 1 
    
    for i in range(numBargroups):
        
        barID = 1
        
        for j in range(barSize):
            
            out.write('>%s-%d-%d\n' % (header, bargroupCount, barID))
            out.write(genomeSeqs[i][j] + '\n')

            barID += 1

        bargroupCount += 1
         

if __name__ == "__main__":
    
    readLen = 150       # Length of read (single-end)
    
    barSize = 200       # Size of bargroup

    genera = {
    'Propionibacterium': 4108,
    'Alteromonas': 7928,
    'Enterobacter': 236,
    'Haemophilus': 119,
    'Neisseria': 77,
    'Staphylococcus': 312,
    'Streptococcus': 754,
    'Bacillus': 1390,
    'Escherichia': 172,
    'Delftia': 274,
    'Pseudomonas': 157,
    'Stenotrophomonas': 187
    }
    
    fastaDir = '/home/ben/ace-seq/ace11.2/virulence_control/genera/'
    
    # Generate FASTA files containing barcode groups for the selected genera:
    
    for genus in genera:
    
        fasta2bargroups(genus + '.fna', readLen, barSize, genera[genus])
    
    
    
    