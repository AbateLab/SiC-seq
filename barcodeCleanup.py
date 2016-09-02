# -*- coding: utf-8 -*-
"""
@author: Ben Demaree

Modified by Ben 5/26/2016
"""

from __future__ import division
import os
import os.path
import collections
import sqlite3
import shlex
import subprocess
import argparse
import timeit
from itertools import izip

####################################################################################  

Qscore = dict((chr(i),i-33) for i in range(33,90))      # Create quality score dictionary

def qualityControl (R1, R2, RI, out1, out2, barQ = 30, readQ = 30, readL = 20, limit = 0, 
                    mask = False, filterComplexity = True, truncate = False, 
                    lowQFracTh = 0.2, readCompTh = 0.85):
# This function takes raw Illumina FASTQ files and filters the reads by quality score. 
    
    R1 = open(R1, 'r')          # Open input FASTQ files
    RI = open(RI, 'r')
    
    out1 = open(out1, 'w')      # Create output FASTQ files
    
    # For paired-end analysis:
    if not SE:
        R2 = open(R2, 'r')
        out2 = open(out2, 'w') 
        
    # Create a counter for the bargroups where keys are barcodes and values are counts:
    bargroups = collections.Counter()       
    
    num = 0                 # Number of barcodes examined
    filteredBar = 0         # Number of barcodes that failed filter
    filteredReads = 0       # Number of reads with excessive low-quality bases (>lowQFracTh)
    lowCReads = 0           # Number of reads with low-complexity (>readCompTh)
    passed = 0              # Number of barcodes that passed filter
    desynch = 0             # Number of desynchronizations between the RI/R1/R2 files

    for line in RI:     # Go through index FASTQ line-by-line
        
        if limit != 0:          # If we set a limit, stop now
            if num >= limit:
                break
        
        num += 1                # Increment read counter
        passFilter = True       # Resets boolean for determining if we should skip this read
    
        firstLine = line.strip()        # Extract first line of index read
        
        if firstLine[0] != "@":         # Ignores line if it doesn't start with an @ character
            print "Index file line ignored because no @"    
            print firstLine
            continue                    # Try the next line
        
        # Look in index read for barcode:
        else:
            ID = firstLine                  # Extract entire cluster ID
            barcode = RI.next().strip()     # Extract the 15 bp barcode
            RI.next()                       # Ignore comment line
            barQual = RI.next().strip()     # Extract quality score
            
            # If one barcode base fails quality check, go on to next read:             
            for base in barQual:
                
                if Qscore[base] < barQ:
                    passFilter = False      # Set boolean to False to mark failed barcode
                    filteredBar += 1
                    break
                
        if not passFilter:      # Discard barcode if it does not pass filter and return to index file
            continue                    
        

        # If barcode passes quality filter, then look at the read to see if cluster IDs match:
        else:
            
            passed += 1     # Increment counter for passed barcodes
            
            for line in R1:     # For the first set of paired-end reads in R1
                
                if line[0] != "@":      # Ignore line if it doesn't start with an @ character
                    print "R1 file line ignored because no @"  
                    print line
                    continue            # Try the next line
                
                elif line[21:44] == ID[21:44]:    # Check that the cluster IDs match
                    # If match, extract read data:                    
                    seq1 = R1.next().strip()    # Extract the sequence
                    R1.next()                   # Ignore comment line
                    qual1 = R1.next().strip()   # Extract quality score
                    
                    # Analyze read quality and mask/truncate (optional):
                    seq1, qual1, lowQFrac1 = maskLowQ(seq1, qual1, readQ, mask, truncate)
                    
                    break       # Advance to reads in R2
                
                else: 
                    desynch += 0.5      # Increment desynch counter if cluster IDs do not match
                    R1.next()
                    R1.next()
                    R1.next()
            
            if not SE:
                
                for line in R2:     # For PE reads, we need to look through R2
                    
                    if line[0] != "@":      # Ignore line if it doesn't start with an @ character
                        print "R2 file line ignored because no @"                     
                        print line                            
                        continue
                    
                    elif line[21:44] == ID[21:44]:  # Check that the cluster IDs match
                        # If match, extract read data:                      
                        seq2 = R2.next().strip()    # Extract the sequence
                        R2.next()                   # Ignore comment line
                        qual2 = R2.next().strip()   # Extract quality score
                        
                        # Analyze read quality and mask/truncate (optional):
                        seq2, qual2, lowQFrac2 = maskLowQ(seq2, qual2, readQ, mask, truncate)
                        break
                    
                    else: 
                        desynch += 0.5      # Increment desynch counter if cluster IDs do not match
                        R2.next()
                        R2.next()
                        R2.next()
                    
            # Throw out pair if read 1 or read 2 are shorter than readL:
            if SE:
                if len(seq1) < readL:                           # SE
                    num += 1
                    filteredReads += 1
                    continue
            else:
                if (len(seq1) < readL) or (len(seq2) < readL):  # PE
                    num += 2
                    filteredReads += 2
                    continue
                                
            # Throw out pair if read 1 or read 2 contain too many low-quality bases:
            if SE:
                if lowQFrac1 > lowQFracTh:                                  # SE
                    num += 1
                    filteredReads += 1
                    continue
            else:
                if (lowQFrac1 > lowQFracTh) or (lowQFrac2 > lowQFracTh):    # PE
                    num += 2
                    filteredReads += 2
                    continue
            
            # Check for sequence complexity:
            if filterComplexity and SE:                 # SE
                seq1C = readComplexity(seq1)

                # Throw out read if it lacks complexity:                        
                if seq1C > readCompTh:
                    num += 1
                    lowCReads += 1
                    continue
                
            if filterComplexity and not SE:             # PE
                seq1C = readComplexity(seq1)
                seq2C = readComplexity(seq2)

                # Throw out pair if read 1 or read 2 lack complexity:                        
                if (seq1C > readCompTh) or (seq2C > readCompTh):
                    num += 2
                    lowCReads += 2
                    continue
                
            # If passing all filters, write the paired reads to new FASTQ files:
            out1.write("@%s-1\n%s\n+\n%s\n" % (barcode, seq1, qual1))
            if not SE:
                out2.write("@%s-2\n%s\n+\n%s\n" % (barcode, seq2, qual2))
            
            # Increment the bargroup counter if it exists in the dictionary:
            if barcode in bargroups:        
                bargroups[barcode] += 1
            
            # If not, make a new entry in the dictonary and start the count at 1:
            else:
                bargroups[barcode] = 1
                            
    print "%d barcodes analyzed" % num    
    print "%d barcodes passed filter" % passed
    print "%d barcodes failed filter" % filteredBar
    print "%d total reads with >%0.1f%% low-quality bases discarded" % (filteredReads, lowQFracTh * 100)
    print "%d total reads with >%0.1f%% identical base calls discarded" % (lowCReads, readCompTh * 100)
    print "%d times desynchronized" % desynch
    print "%d barcode groups identified" % len(bargroups)
    
    R1.close()      # Close raw FASTQ input files            
    RI.close()
    
    out1.close()    # Close output files
    
    if not SE:
        R2.close()
        out2.close()
    
    return bargroups

def maskLowQ (seq, qual, readQ, mask = True, truncate = False):
# This function finds bases in reads with a Qscore < readQ and replaces them 
# with a placeholder 'N'
# If truncate == True, the read will be cut at the first low-quality base
    
    # Find indices of low-quality bases:
    lowQbases = [i for i in range(len(qual)) if Qscore[qual[i]] < readQ]    
    
    lowQFrac = len(lowQbases) / len(seq)      # Fraction of low-quality bases

    if mask:            # With masking
        
        for i in lowQbases:
            
            seq = seq[:i] + 'N' + seq[i+1:]         # Replace low-quality bases with 'N'
                 
    elif truncate:      # With truncation
        
        if len(lowQbases) != 0:
            seq = seq[:min(lowQbases)]      # Truncate the sequence
            qual = qual[:min(lowQbases)]    # Truncate the quality

    return seq, qual, lowQFrac
    
def readComplexity (seq):
# Returns the fraction of base calls in a read corresponding to a single base.

    counts = [seq.count('A'), seq.count('T'), seq.count('C'), seq.count('G'), seq.count('N')]
    
    complexity = max(counts) / len(seq)         # Calculate complexity
    
    return complexity

def findOrphans (bargroups, orphans, barCutoff = 100, save = False):
# Finds barcodes that have less than 'cutoff' reads as orphans and stores them 
# in a txt file if save == True.
    
    bargroupsCutoff = {}    # Create new counter for bargroups passing cutoff
    
    failedCount = 0         # Count number of bargroups failing cutoff
    passedCount = 0         # Count number of bargroups passing cutoff
    
    if save:    # With save option

        orphans = open("orphans.txt", 'w')     # Create output orphans.txt file
        
        for barcode in bargroups:       # Look through all bargroups
            
            if bargroups[barcode] < barCutoff:      # Bargroups failing cutoff   
                # Write failed barcodes to file along with bargroup size:
                orphans.write(barcode + "-" + str(bargroups[barcode]) + "\n")
                failedCount += 1                    # Increment failed bargroup count
                
            else:       # Bargroups passing cutoff
                # Put passing bargroups in new dictionary with values from old bargroups counter:
                bargroupsCutoff[barcode] = "%d-%d" % (bargroups[barcode], bargroups[barcode])
                passedCount += 1    # Increment passed bargroup count           
            
        print "%d bargroups passed cutoff threshold of %d" % (passedCount, barCutoff)
        print "%d bargroups failed cutoff and were saved to orphans.txt" % failedCount
        
        orphans.close()
        
    else:       # Without save option
        
        for barcode in bargroups:       # Look through all bargroups
            
            if bargroups[barcode] >= barCutoff:     # Bargroups passing cutoff
                # Put bargroup in new dictionary with values from old bargroups counter:
                bargroupsCutoff[barcode] = "%d-%d" % (bargroups[barcode], bargroups[barcode])
                passedCount += 1    # Increment passed bargroup count
                
            else:       # Bargroups failing cutoff
                failedCount += 1    # Increment failed bargroup count 
                
        print "%d bargroups passed cutoff threshold of %d reads." % (passedCount, barCutoff)
        print "%d bargroups failed cutoff and were discarded." % failedCount
    
    return bargroupsCutoff

def bargroupID (bargroups, in1, in2, out1, out2):
# This function appends the bargroup ID and bargroup size to each read and 
# saves the bargroups passing cutoff filter to new FASTQ files.

    # Create output FASTQ files:
    out1 = open(out1, 'w')
    if not SE:
        out2 = open(out2, 'w')
    
    # Open input FASTQ files:
    in1 = open(in1, 'r')
    if not SE:
        in2 = open(in2, 'r')
    
    readCount = 0       # Count of reads written to file
    
    for line in in1:    # Look through the input file line-by-line
        
        barcode = line.strip()[1:16]    # Extract 15 bp barcode from R1 input file
        if not SE:
            in2.next()                  # Skip the ID line in R2 input
        
        if barcode in bargroups:        # If reads in FASTQ belong to a bargroup passing filter
            
            ID_Size = bargroups[barcode].split("-")     # Split ID line
            barID = int(ID_Size[0])                     # Extract barcode ID
            barSize = int(ID_Size[1])                   # Extract bargroup size
            
            # For R1 FASTQ:
            seq1 = in1.next().strip()       # Extract sequence 1 
            in1.next()                      # Skip comment line
            qual1 = in1.next().strip()      # Extract quality 1
            
            out1.write("@%s-1-%d-%d\n%s\n+\n%s\n" % (barcode, barID, barSize, seq1, qual1))     # Write to file
            
            if not SE:
                # For R2 FASTQ:
                seq2 = in2.next().strip()       # Extract sequence 2
                in2.next()                      # Skip comment line
                qual2 = in2.next().strip()      # Extract quality 2
                
                out2.write("@%s-2-%d-%d\n%s\n+\n%s\n" % (barcode, barID, barSize, seq2, qual2))     # Write to file
            
            barID -= 1      # Decrease barcode ID after saving read
                
            # Save new barcode ID and bargroup size in dictionary:
            bargroups[barcode] = "%d-%d" % (barID, barSize)
                
            readCount += 1      # Increment read counter
                
            if (readCount % 1e6 == 0):
                print "%d reads ID'ed so far..." % readCount
            
        # If read belongs to a bargroup not passing filter, skip ahead 3 lines 
        # to the next barcode ID in FASTQ files: 
        else:
            try:    
                in1.next(); in1.next(); in1.next()
                if not SE:
                    in2.next(); in2.next(); in2.next()
            
            except (StopIteration):     # When the end of the file is reached
                print "%d reads tagged with ID." % readCount
                break
    
    # Sanity check to make sure all barcodes are accounted for:
    for barcode in bargroups:
        
        if barID != 0:
            print "Barcode counter must equal 0 for all bargroups!"
            raise SystemExit    
    
    in1.close()     # Close input FASTQ files
    if not SE:
        in2.close()
    
    out1.close()    # Close output FASTQ files
    if not SE:
        out2.close()
    
def sendToBowtie(samfile, bt2index, FASTQ1, FASTQ2):
# Sends the bargroups passing cutoff to bowtie2 for alignment.
    
    if not SE:
        bowtieCmd = "bowtie2 -x %s -1 %s -2 %s -S %s --no-sq -p 5 --un-conc unaligned_R%%.fastq" \
        % (bt2index, FASTQ1, FASTQ2, samfile)
        
    else:
        bowtieCmd = "bowtie2 -x %s -U %s -S %s --no-sq -p 5 --un-conc unaligned_R1.fastq" \
        % (bt2index, FASTQ1, samfile)
    
    align = subprocess.call(shlex.split(bowtieCmd))
    
    print "SAM file %s saved!" % samfile
        
def toSQLiteFASTQ (FASTQ1, FASTQ2, sqliteDB):
# This function saves the cleaned-up bargroups in FASTQ files to an SQLite database.

    FASTQ1 = open(FASTQ1, 'r')          # Open input FASTQ files
    if not SE:
        FASTQ2 = open(FASTQ2, 'r')
    
    readCount = 0       # Create counter for reads handled
    
    # Initialize the SQLite DB:
    sql = sqlite3.connect(sqliteDB, isolation_level='Exclusive')
    s = sql.cursor()            # Create cursor
    # These settings allow DB to be stored in memory and written periodically to the disk: 
    s.execute('PRAGMA synchronous = 0')
    s.execute('PRAGMA journal_mode = OFF')
    
    lineCount = 0
    
    if not SE:
        # Create a table in the DB for storing the reads:
        s.execute('''CREATE TABLE reads
                    (barcode TEXT, barID INT, barSize INT, seq1 TEXT, seq2 TEXT, 
                     startPos1 INT, startPos2 INT, species TEXT)''')
                
        for line1, line2 in izip(FASTQ1, FASTQ2):     # Iterate through lines of the FASTQ files together
    
            lineCount += 1
            
            if lineCount % 4 == 1:
            
                # FASTQ for Read 1:
                IDLine = line1.split("-")       # Split up the ID line
                barcode = IDLine[0][1:16]       # Extract the 15 bp barcode
                barID = int(IDLine[2])          # Extract the barcode ID
                barSize = int(IDLine[3])        # Extract the bargroup size
                
                # FASTQ for Read 2:
                IDLine2 = line2.split("-")      # First line in FASTQ 2 (identical to FASTQ 1)
                barcode2 = IDLine2[0][1:16]     # Extract the 15 bp barcode
                
                # Check that the barcodes match:
                if barcode2 != barcode:
                    print 'Mismatch between input FASTQ files during export to DB!'
                    raise SystemExit
                
            elif lineCount % 4 == 2:
                seq1 = line1.strip()     # Extract sequence 1
                seq2 = line2.strip()     # Extract sequence 2
    
            elif lineCount % 4 == 3:
                continue                # Skip comment line
                
            else:
                
                # Create a row in the DB for the extracted read information
                # By default, we set the species name to 'unal' and update this
                # information after species identification:
                s.execute('''INSERT INTO reads (barcode, barID, barSize, 
                seq1, seq2, species, startPos1, startPos2) values (?, ?, ?, ?, ?, ?, ?, ?)''',
                (barcode, barID, barSize, seq1, seq2, "unal", -1, -1))
                
                readCount += 1      # Increment read counter
                    
                # We periodically need to flush the memory and write DB entries to disk:
                if (readCount % 1e6 == 0):
                    sql.commit()        # Write reads to persistent DB
                    print "%d reads exported to database so far..." % readCount
            
        print "%d reads exported to database!" % readCount
    
        FASTQ1.close()      # Close input files
        FASTQ2.close()   
                    
    else:
        # Create a table in the DB for storing the reads:
        s.execute('''CREATE TABLE reads
                    (barcode TEXT, barID INT, barSize INT, seq1 TEXT, 
                     startPos1 INT, species TEXT)''')
                
        for line1 in FASTQ1:     # Iterate through lines of the FASTQ file
    
            lineCount += 1
            
            if lineCount % 4 == 1:
            
                # FASTQ for Read 1:
                IDLine = line1.split("-")       # Split up the ID line
                barcode = IDLine[0][1:16]       # Extract the 15 bp barcode
                barID = int(IDLine[2])          # Extract the barcode ID
                barSize = int(IDLine[3])        # Extract the bargroup size
                
            elif lineCount % 4 == 2:
                seq1 = line1.strip()     # Extract sequence 1
                
            elif lineCount % 4 == 3:
                continue                # Skip comment line
                
            else:
                
                # Create a row in the DB for the extracted read information
                # By default, we set the species name to 'unal' and update this
                # information after species identification:
                s.execute('''INSERT INTO reads (barcode, barID, barSize, 
                seq1, species, startPos1) values (?, ?, ?, ?, ?, ?)''',
                (barcode, barID, barSize, seq1, "unal", -1))
                
                readCount += 1      # Increment read counter
                    
                # We periodically need to flush the memory and write DB entries to disk:
                if (readCount % 1e6 == 0):
                    sql.commit()        # Write reads to persistent DB
                    print "%d reads exported to database so far..." % readCount
            
        print "%d reads exported to database!" % readCount
    
        FASTQ1.close()      # Close input file

    sql.commit()        # Commit db changes to disk one last time
    sql.close()         # Close the db
    
    print "SQLite database saved!"
    
def toSQLiteSAM (samfile, sqliteDB, organisms):
# This function exports the bowtie2-aligned sequences in SAM format to an SQLite database.

    sam = open(samfile, 'r')          # Open input FASTQ files
    
    readCount = 0       # Create counter for reads handled
     
    # Initialize the SQLite DB:
    sql = sqlite3.connect(sqliteDB, isolation_level='Exclusive')
    s = sql.cursor()            # Create cursor
    # These settings allow DB to be stored in memory and written periodically to the disk: 
    s.execute('PRAGMA synchronous = 0')
    s.execute('PRAGMA journal_mode = OFF')

    
    for organism in organisms:
        exec("%sCount = 0" % organism)

    unalignedCount = 0      # Number of reads not aligned to a reference
    alignedCount = 0        # Number of reads aligned to a reference
    mismatchCount = 0       # Number of reads where species do not match
    
    # For paired-end reads:
    if not SE:    
        # Create a table in the DB for storing the reads:
        s.execute('''CREATE TABLE reads
                    (barcode TEXT, barID INT, barSize INT, seq1 TEXT, seq2 TEXT, 
                    startPos1 INT, startPos2 INT, species TEXT)''')
        for line in sam:
                
            if line[0] == "@":  # Ignore header lines
                continue
            
            # Split each SAM file line for parsing:
            splitLine1 = line.split("\t")
            
            idSplit1 = splitLine1[0].split("-")        
            barcode1 = idSplit1[0]      # Barcode 1
            pairNo1 = idSplit1[1]       # Read-pair number (1 or 2)
            barID = idSplit1[2]         # Read-pair ID (1, 2,...barSize)
            barSize = idSplit1[3]       # Number of read-pairs in bargroup
            
            seq1 = splitLine1[9]                # Sequence 1
            ref1 = splitLine1[2]                # Reference seq1 aligns to
            flag1 = int(splitLine1[1])          # Bitwise flag 1
            unmapped1 = int(bin(flag1)[-3])     # Mapping bit (1: not mapped; 0: mapped)
            
            # Go to next line (should be next mate in read-pair):
            splitLine2 = sam.next().split("\t")
            barcode2 = splitLine2[0].split("-")[0]      # Barcode 2
            
            # Make sure barcodes match:
            if barcode1 != barcode2:
                print 'Barcodes do not match in SAM file! Check inputs.'
                raise SystemExit
                
            seq2 = splitLine2[9]                # Sequence 2
            ref2 = splitLine2[2]                # Reference seq2 aligns to
            flag2 = int(splitLine2[1])          # Bitwise flag 2
            unmapped2 = int(bin(flag2)[-3])     # Mapping bit (1: not mapped; 0: mapped)     
            
            # Check FLAG field in SAM file for mapping status:
            
            if unmapped1 == 1 and unmapped2 == 1:       # Both reads unaligned
                species = "unal"
                startPos1 = -1
                startPos2 = -1
                unalignedCount += 2
            
            elif unmapped1 == 0 and unmapped2 == 1:     # Only R1 aligned
                org1 = ref1.split("-")[0]
                species = org1
                startPos1 = int(splitLine1[3])
                startPos2 = -1
                
            elif unmapped1 == 1 and unmapped2 == 0:     # Only R2 aligned
                org2 = ref2.split("-")[0]
                species = org2
                startPos1 = -1
                startPos2 = int(splitLine2[3])
            
            else:       # Both reads aligned
                org1 = ref1.split("-")[0]
                org2 = ref2.split("-")[0]
                
                # Detecting mismatches:
                if org1 != org2:        # Both reads align, but to different organisms
                    mismatchCount += 2  
                    continue            # Don't add mismatches to the database
                
                # When both mates align to the same reference:
                else:
                    species = org1
                    startPos1 = int(splitLine1[3])
                    startPos2 = int(splitLine2[3])
                    alignedCount += 2
                    exec("%sCount += 1" % org1)
                    
            # Add the paired reads to the table:
            
            if pairNo1 == 1:    # Ensure mates are added to the correct field in the db
                
                s.execute('''INSERT INTO reads (barcode, barID, barSize, 
                seq1, seq2, species, startPos1, startPos2) values (?, ?, ?, ?, ?, ?, ?, ?)''',
                (barcode1, barID, barSize, seq1, seq2, species, startPos1, startPos2))
                
            else:
                
                s.execute('''INSERT INTO reads (barcode, barID, barSize, 
                seq1, seq2, species, startPos1, startPos2) values (?, ?, ?, ?, ?, ?, ?, ?)''',
                (barcode1, barID, barSize, seq2, seq1, species, startPos2, startPos1))
                
            readCount += 2      # Increment read counter
                
            # We periodically need to flush the memory and write DB entries to disk:
            if (readCount % 2e6 == 0):
                sql.commit()        # Write reads to persistent DB
                print "%d total reads exported to database so far..." % readCount
    
    # For single-end reads:          
    else:    
        # Create a table in the DB for storing the reads:
        s.execute('''CREATE TABLE reads
                    (barcode TEXT, barID INT, barSize INT, seq1 TEXT, 
                    startPos1 INT, species TEXT)''')
        for line in sam:
                
            if line[0] == "@":  # Ignore header lines
                continue
            
            # Split each SAM file line for parsing:
            splitLine1 = line.split("\t")
            
            idSplit1 = splitLine1[0].split("-")        
            barcode1 = idSplit1[0]      # Barcode 1
            barID = idSplit1[2]         # Read-pair ID (1, 2,...barSize)
            barSize = idSplit1[3]       # Number of read-pairs in bargroup
            
            seq1 = splitLine1[9]                # Sequence 1
            ref1 = splitLine1[2]                # Reference seq1 aligns to
            flag1 = int(splitLine1[1])          # Bitwise flag 1
            unmapped1 = int(bin(flag1)[-3])     # Mapping bit (1: not mapped; 0: mapped)
            
            # Check FLAG field in SAM file for mapping status:
            
            if unmapped1 == 1:      # Read is unaligned
                species = "unal"
                startPos1 = -1
                unalignedCount += 1
            
            elif unmapped1 == 0:    # Read is aligned
                org1 = ref1.split("-")[0]
                species = org1
                startPos1 = int(splitLine1[3])
                alignedCount += 1
                exec("%sCount += 1" % org1)
                    
            # Add the single-end read to the table:
            s.execute('''INSERT INTO reads (barcode, barID, barSize, 
            seq1, species, startPos1) values (?, ?, ?, ?, ?, ?)''',
            (barcode1, barID, barSize, seq1, species, startPos1))
                
            readCount += 1      # Increment read counter
                
            # We periodically need to flush the memory and write DB entries to disk:
            if (readCount % 2e6 == 0):
                sql.commit()        # Write reads to persistent DB
                print "%d total reads exported to database so far..." % readCount
        
    print "%d reads exported to database!" % readCount

    sam.close()         # Close SAM file
    
    sql.commit()        # Commit db changes to disk one last time
    sql.close()         # Close the db
    
    print "SQLite database saved!"
    
    totalUnique = unalignedCount    # Total unique reads (add in alignments next)
    
    for organism in organisms:
        exec("totalUnique += %sCount" % organism)  
    
    print "%d total reads examined. Of these..." % readCount
    
    for organism in organisms:
        
        exec("percentAligned = %sCount / totalUnique * 100" % organism)
        
        exec('''print("\t" + str(%sCount) + " (" + str(round(percentAligned, 2)) + 
        "%%) aligned to the %s genome")''' % (organism, organism))    
    
    print "\t%d (%0.2f%%) did not align\n" % (unalignedCount, unalignedCount / totalUnique * 100)
    
    print "There were %d mismatched reads" % mismatchCount    
    
if __name__ == "__main__":
 
####################################################################################  
# Step 0: Parse input arguments and set script variables
#################################################################################### 
 
    ### Setup argument parser ###
     
    parser = argparse.ArgumentParser(description='''
    
    Clean up FASTQ files containing sequencing information for an ACE-seq 
    experiment. This script filters reads by basecall quality and barcode group
    size and exports the results to an SQLite database. Runs on paired-end
    FASTQ files by default (for SE, specify --SE flag).
    
    Before running, set the 'organisms' list in barcodeCleanup.py to match the
    headers in the FASTA files used to build the bowtie2 index. Use the script
    refLabels.py to relabel the headers of FASTA files.
    
    ***Dependencies***
    sqlite3
    
    ''', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('label', type=str,
                        help='experiment label for this sample')                
                            
    parser.add_argument('-R1', type=str, metavar='FASTQ_file', default='R1.fastq',
                        help='FASTQ file containing Read 1 [default: R1.fastq]') 
                        
    parser.add_argument('-R2', type=str, metavar='FASTQ_file', default='R2.fastq',
                        help='FASTQ file containing Read 2 [default: R2.fastq]')
                        
    parser.add_argument('-RI', type=str, metavar='FASTQ_file', default='RI.fastq',
                        help='FASTQ file containing index reads [default: RI.fastq]')
    
    parser.add_argument('-align', choices=['Y', 'N'], metavar='Y/N', default='Y',
                        help='option to align reads to one or more reference genomes [default: Y]')
                            
    parser.add_argument('-bt2', type=str, metavar='bt2index', default='none',
                        help='bowtie2 index base name [default: experiment label]')
        
    parser.add_argument('-barC', type=int, metavar='int', default=50,
                        help='bargroup cutoff size [default: 50 reads (SE)/read-pairs (PE)]') 
                        
    parser.add_argument('-barQ', type=int, metavar='int', default=20,
                        help='minimum barcode base quality [default: 20]')
                        
    parser.add_argument('-readQ', type=int, metavar='int', default=20,
                        help='minimum read base quality [default: 20]')
                        
    parser.add_argument('-readL', type=int, metavar='int', default=20,
                        help='minimum read length [default: 20]')
    
    parser.add_argument('-lowF', type=float, metavar='frac', default=0.20,
                        help='maximum fraction of low-quality (<readQ) bases in a read [default: 0.20]')
    
    parser.add_argument('-singF', type=float, metavar='frac', default=0.85,
                        help='maximum fraction of a single base type in a read [default: 0.85]')
                        
    parser.add_argument('-complex', choices=['Y', 'N'], metavar='Y/N', default='Y',
                        help='option to filter out low-complexity (>singF) reads [default: Y]')
                        
    parser.add_argument('-truncate', choices=['Y', 'N'], metavar='Y/N', default='N',
                        help='option to truncate read at first low-quality (<readQ) base [default: N]')
                        
    parser.add_argument('-mask', choices=['Y', 'N'], metavar='Y/N', default='N',
                        help='option to mask low-quality (<readQ) bases [default: N]')
                                                
    parser.add_argument('-limit', type=int, metavar='int', default=0,
                        help='maximum number of read-pairs to analyze [default: 0 - no limit]')
                        
    parser.add_argument('--SE', action='store_true', default=False,
                        help='flag for a single-end read analysis') 
    
    args = parser.parse_args()      # Parse arguments
    
    ### Set variables and filenames based on command line inputs ###

    expName = args.label    # Experiment name

    barCutoff = args.barC   # Minimum number of read-pairs in a bargroup
    
    # Raw FASTQ files from the sequencer:
    R1 = args.R1
    R2 = args.R2
    RI = args.RI
    
    # Options for read filtering:    
    barQ = args.barQ            # Minimum barcode quality score
    readQ = args.readQ          # Minimum base quality score in reads
    readL = args.readL          # Minimum read length
    lowQFracTh = args.lowF      # Maximum fraction of low-quality bases in a read
    readCompTh = args.singF     # Maximum abundance of a single base type in a read
    SE = args.SE                # Option to analyze a single-end read
    
    YN = {'Y': True, 'N': False}    # For mapping input arguments to booleans  
    
    align = YN[args.align]                  # Option to align reads to reference
    mask = YN[args.mask]                    # Option to mask low-quality read bases (<readQ)
    truncate = YN[args.truncate]            # Option to truncate read at first low-quality base (<readQ)
    filterComplexity = YN[args.complex]     # Option to filter reads by complexity (>readCompTh)
    
    limit = args.limit           # Maximum number of reads from raw input to parse (0: no limit)
    
    # Base name of pre-built bowtie2 index:
    if args.bt2 == 'none':
        bt2index = args.label
    else:
        bt2index = args.bt2
        
    # Make sure that bowtie2 index exists:
    if not os.path.isfile(bt2index + '.1.bt2') and align:
        print 'Specified bt2 index not found! Please check the file name and path.'
        raise SystemExit
    
    # Make sure DB with same name does not exist:
    sqliteDB = "%s-ACE-%dbarCutoff.db" % (expName, barCutoff)       # File name of sqlite DB
        
    if os.path.isfile(sqliteDB):        # Don't overwrite an existing DB
        print "The specified DB file already exists. Either rename or delete the existing file."
        raise SystemExit

####################################################################################  
# Step 1: filter out low-quality barcodes and reads
####################################################################################   

    # Output intermediate FASTQ files with filtered barcodes.
    # Also create 'bargroups' dictionary containing all bargroup sizes:

    print "Beginning barcode cleanup routine..."
    tic = timeit.default_timer()    # Record start time 
    
    print "Filtering barcode groups by index and read quality score..."
    
    # File names of temporary output FASTQ files (to be deleted):
    out1 = "R1_QC.fastq"
    out2 = "R2_QC.fastq"

    bargroups = qualityControl (R1, R2, RI, out1, out2, barQ, readQ, readL, limit, 
                                mask, filterComplexity, truncate, 
                                lowQFracTh, readCompTh)
    
    # Text file name for storing bargroups failing cutoff:
    orphans = "orphans-%s-%dbarCutoff.txt" % (expName, barCutoff)
    
    # Output dictionary with filtered bargroups (# reads > barCutoff):
    
    print "Finding and removing bargroups with less than %d reads..." % barCutoff
    
    bargroupsCutoff = findOrphans (bargroups, orphans, barCutoff, True)

####################################################################################    
# Step 2: append barcode groups with size information and remove small bargroups
####################################################################################
    
    # File names of FASTQ files prepared for alignment and post-processing:
    FASTQ1 = "%dbarCutoff_%s" % (barCutoff, out1)
    FASTQ2 = "%dbarCutoff_%s" % (barCutoff, out2)    
    
    # Append barcode ID and bargroup size to reads and export to new FASTQ files:
    
    print "Appending barcode ID and bargroup size to reads..."    
    
    bargroupID (bargroupsCutoff, out1, out2, FASTQ1, FASTQ2)
    
    # Delete temporary FASTQ files:
    os.remove(out1)
    if not SE:
        os.remove(out2)

####################################################################################    
# Step 3: align to reference genomes using bowtie2 [if selected]
####################################################################################
    
    if align:
        
        # Known organisms in sample (must match headers in FASTA files):
        organisms = ['bs168', 'yeast', 'se']   
        
        samfile = "%s-%dcutoff.sam" % (expName, barCutoff)  # Name of output SAM file
        
        print "Beginning alignment using bowtie2..."
        
        # Align using bowtie2 and save to SAM file:
        sendToBowtie(samfile, bt2index, FASTQ1, FASTQ2)
    
####################################################################################
# Step 4: add reads to SQLite database
####################################################################################
    
    if align:    
    
        # Export reads from SAM files into an SQLite database stored on disk:
        
        print "Beginning export to SQLite database from SAM file..."    
        
        toSQLiteSAM(samfile, sqliteDB, organisms)
        
    else:
        
        # Export reads from FASTQ files into an SQLite database stored on disk:
        
        print "Beginning export to SQLite database from FASTQ files..."    
        
        toSQLiteFASTQ(FASTQ1, FASTQ2, sqliteDB)
    
    toc = timeit.default_timer()    # Record stop time
    print "%d seconds total processing time." % (toc-tic)
    