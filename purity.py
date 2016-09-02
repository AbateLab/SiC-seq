# -*- coding: utf-8 -*-
"""
Created on Thu May 19 11:39:58 2016

@author: Ben Demaree
"""

### Misc. imports
from __future__ import division
import sqlite3
from scipy.stats import mode
import argparse

### Pandas imports
import pandas as pd
from sqlalchemy import create_engine

####################################################################################

def computePurity (table):
# Calculates the purity of the bargroups in a specified database table.
    
    df = pd.read_sql_query('select barcode,species from "%s" where species!="unal" order by barcode' % table, db)
    
    barcode = list(df['barcode'])
    species = list(df['species'])
    
    print '%d read-pairs are aligned to a reference.' % len(barcode)
    print '%d barcode groups contain at least one aligned read-pair.' % len(set(barcode))
    print 'Calculating barcode group purities from aligned reads...'
       
    purities = []
    sizes = []
    mostCommonSpecies = []
    barcodes = []
    
    last = 0
    
    for i in range(len(barcode) - 1):
        
        # Continue while still in the same bargroup:
        if barcode[i] == barcode[i + 1]:
            continue
        
        # Once we reach the end of the bargroup:
        else:
            barSize = i - last + 1      # Number of aligned reads in the group   
            
            sizes += [barSize]          # Add to sizes list
            
            barcodes += [barcode[i]]    # Add to barcodes list
            
            modeCount = mode(species[last:(i + 1)])     # Find mode of species (i.e. most common)

            mostCommonSpecies += [modeCount[0][0]]      # Add to most common species list
        
            purities += [modeCount[1][0] / barSize]     # Divide count of most common by count of aligned
        
            last = i + 1    # Set bargroup size counter to current position
    
    # Add last bargroup to the list:
    
    barSize = len(barcode) - last         
    
    sizes += [barSize]
    
    barcodes += [barcode[-1]]

    modeCount = mode(species[last:len(barcode)])

    mostCommonSpecies += [modeCount[0][0]]

    purities += [modeCount[1][0] / barSize]
        
    print 'Barcode group purities calculated!'
        
    return purities, sizes, mostCommonSpecies, barcodes
    
def insertPurities (purities, sizes, mostCommonSpecies, barcodes, table):
# Inserts purity and species information into the sqlite database.

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
    
    # Add columns for the purity and most common species (if not already done):
    try:    # purity
        s.execute('alter table "%s" add column purity real' % table)        
    except sqlite3.Error:
        pass
    
    try:    #mostCommonSpecies
        s.execute('alter table "%s" add column mostCommonSpecies text' % table)        
    except sqlite3.Error:
        pass
            
    # Insert purity and most common species information:
    for i in range(len(barcodes)):

        # If there are not enough aligned reads to assign a purity score:
        if sizes[i] < minReads:
            purities[i] = -1
            mostCommonSpecies[i] = ''
            
        s.execute('update "%s" set purity="%s",mostCommonSpecies="%s" where barcode="%s"' %
        (table, purities[i], mostCommonSpecies[i], barcodes[i]))
        
        if (i % 1e4 == 0 and i != 0):
            sql.commit()        # Write reads to persistent DB
            print "%d purity scores exported to database so far..." % i
                
    sql.commit()        # Commit db changes to disk one last time
    sql.close()         # Close the db
    
    print '%d barcode group purities inserted into database!' % len(barcodes)

if __name__ == "__main__":
    
    ### Setup argument parser ###
     
    parser = argparse.ArgumentParser(description="""
    Calculates the purity within barcode groups using alignment data. Purity 
    scores are added to the specified SQLite database.
    
    A minimum number of aligned reads can be required to assign a purity score
    to a bargroup. Groups that do not meet the threshold are assigned a
    purity score of -1.
    
    ***Dependencies***
    sqlite3
    
    @author: Ben Demaree
    """, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('db', metavar='<db_filename>', type=str,
                        help='SQLite input database filename')
                        
    parser.add_argument('-m', type=int, metavar='min_reads', default=0,
                        help='minimum number of aligned reads/read-pairs per bargroup [default: 0 - no minimum]')
                        
    parser.add_argument('-t', metavar='table', type=str, default='reads',
                        help='table in database to analyze [default: reads]')
        
    args = parser.parse_args()      # Parse arguments
        
    sqliteDB = args.db              # DB filename
    table = args.t                  # Table in DB to analyze
    minReads = args.m               # Minimum number of reads in a bargroup
    
    expName = sqliteDB.split('.')[0]    # Universal filename label
    
    db = create_engine('sqlite:///' + sqliteDB)     # Create DB interface
    
    ### Calculate purities and insert into database ###
    
    purities, sizes, mostCommonSpecies, barcodes = computePurity(table)
    
    insertPurities(purities, sizes, mostCommonSpecies, barcodes, table)



        
        
    
    
    
    

    