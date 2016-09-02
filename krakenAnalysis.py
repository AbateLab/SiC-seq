# -*- coding: utf-8 -*-
"""
Created on Wed Jun 01 11:56:57 2016

@author: Ben

"""
############## Modules ##############

from __future__ import division
import sqlite3
import shlex
import subprocess
import numpy as np
from collections import Counter
import math
import argparse
from scipy.stats import mode

### Matplotlib imports
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.family'] = 'Ubuntu'

### Pandas imports
import pandas as pd
from sqlalchemy import create_engine

############## Functions ##############

def getNumClassified (table):
# Inserts the number of classified reads for each bargroup into the database.
    
    # Pull all data from reads with phylogenetic data in this column:
    df = pd.read_sql_query('select barcode,classified from "%s" order by barcode' % table, db)
    
    barcode = list(df['barcode'])       # Barcodes from classified reads
    classified = list(df['classified']) # Classification status of a read
    
    numClassified = []          # Number of classified reads in a bargroup
    barcodes = []               # Barcodes analyzed
    
    last = 0    # Index for beginning for current group
    
    for i in range(len(barcode) - 1):
        
        # Continue while still in the same bargroup:
        if barcode[i] == barcode[i + 1]:
            continue
        
        # Once we reach the end of the bargroup:
        else:
            barSize = i - last + 1      # Number of  reads in the group (this layer)
            
            barcodes += [barcode[i]]    # Add to barcodes list
            
            numClassified += [classified[last:(i + 1)].count('Y')]    # Number of classified reads
        
            last = i + 1    # Set bargroup size counter to current position
    
    # Add last bargroup to the list:
    
    barSize = len(barcode) - last
    
    barcodes += [barcode[-1]]
    
    numClassified += [classified[last:len(barcode)].count('Y')]
    
    # Import into db:
    
    sql = sqlite3.connect(sqliteDB, isolation_level = 'Exclusive')
    s = sql.cursor()    # Create cursor
    
    # These settings allow DB to be stored in memory and written periodically to the disk: 
    s.execute('PRAGMA synchronous = 0')
    s.execute('PRAGMA journal_mode = OFF')
    
    # Add column for numClassified:
    try:
        s.execute('alter table "%s" add column numClassified int' % table)
        
    except sqlite3.Error:
        pass
    
    # Insert into db:
    for i in range(len(barcodes)):
        s.execute('update "%s" set numCLassified="%s" where barcode="%s"' %
        (table, numClassified[i], barcodes[i]))
        
        if (i % 1e3 == 0 and i != 0):
            sql.commit()        # Write reads to persistent DB
            print "Statistics for %d barcode groups exported to database so far..." % i
            
    sql.commit()        # Commit db changes to disk
    sql.close()         # Close the db

def computeStats (table, column, classifiedSizes):
# Calculates the statistics for Kraken data in a specified table.
    
    # Pull all data from reads with phylogenetic data in this column:
    df = pd.read_sql_query('select barcode,barSize,numClassified,"%s" from "%s" where "%s"!="" order by barcode' % (column, table, column), db)
    
    barcode = list(df['barcode'])       # Barcodes from classified reads
    sizeAll = list(df['barSize'])       # Size of ENTIRE group
    columndb = list(df[column])         # Column of interest
    numClassified = list(df['numClassified'])   # Number of classified reads in bargroup
       
    purities = []               # Purity for given column
    fracClassifiedLayer = []    # Fraction of classified reads for this phylogenetic layer
    fracClassified = []         # Fraction of classified reads in a bargroup
    mostCommon = []             # Most common entry in given column
    barcodes = []               # Barcodes analyzed
    
    last = 0    # Index for beginning for current group
    
    for i in range(len(barcode) - 1):
        
        # Continue while still in the same bargroup:
        if barcode[i] == barcode[i + 1]:
            continue
        
        # Once we reach the end of the bargroup:
        else:
            barSize = i - last + 1      # Number of classified reads in the group (this layer)

            fracClassifiedLayer += [barSize / numClassified[i]]   # % classified (specified layer only)
            
            fracClassified += [barSize / sizeAll[i]]    # % classified (all reads)
            
            barcodes += [barcode[i]]    # Add to barcodes list
            
            modeCount = mode(columndb[last:(i + 1)])    # Find mode of column (i.e. most common)

            mostCommon += [modeCount[0][0]]             # Add to most common species list
        
            purities += [modeCount[1][0] / barSize]     # Divide count of most common by count of aligned
        
            last = i + 1    # Set bargroup size counter to current position
    
    # Add last bargroup to the list:
    
    barSize = len(barcode) - last
    
    fracClassifiedLayer += [barSize / numClassified[-1]]
    
    fracClassified += [barSize / sizeAll[-1]]
    
    barcodes += [barcode[-1]]

    modeCount = mode(columndb[last:len(barcode)])

    mostCommon += [modeCount[0][0]]

    purities += [modeCount[1][0] / barSize]
        
    print 'Data computed for column %s!' % column
        
    return purities, fracClassified, fracClassifiedLayer, mostCommon, barcodes
    
def insertStats (table, column, purities, fracClassified, fracClassifiedLayer, mostCommon, barcodes):
# Creates new columns in the table with the computed statistics
    
    print 'Inserting phylogenetic statistics into database...'
    
    sql = sqlite3.connect(sqliteDB, isolation_level = 'Exclusive')
    s = sql.cursor()    # Create cursor
    
    # These settings allow DB to be stored in memory and written periodically to the disk: 
    s.execute('PRAGMA synchronous = 0')
    s.execute('PRAGMA journal_mode = OFF')
    
    # Classified reads fraction (all layers) - only need to do this once:
    if column == 'domain':
        try:
            s.execute('alter table "%s" add column fracClassified float default 0' % table)    
        except sqlite3.Error:
            pass
        
    # Add index for barcodes (if not already done):
    try:   
        # Create index on barcodes:
        s.execute('create index indexMain on "%s" (barcode)' % table)
    
    except sqlite3.Error:
        pass
    
    # Add layer-specific columns:
    try:
        s.execute('alter table "%s" add column "%s" float' % (table, 'purity' + column))
        s.execute('alter table "%s" add column "%s" float' % (table, 'fracClassified' + column))
        s.execute('alter table "%s" add column "%s" text' % (table, 'mostCommon' + column))
        
    except sqlite3.Error:
        pass
      
    # Insert phylogenetic statistics data:
    for i in range(len(barcodes)):
        
        # Insert bargroup classified fraction with domain data only:
        if column == 'domain':            
            s.execute('update "%s" set fracClassified="%s","%s"="%s","%s"="%s","%s"="%s" where barcode="%s"' %
            (table, fracClassified[i], 'purity' + column, purities[i], 'fracClassified' + column, fracClassifiedLayer[i], 'mostCommon' + column, mostCommon[i], barcodes[i]))
            
            if (i % 1e3 == 0 and i != 0):
                sql.commit()        # Write reads to persistent DB
                print "Statistics for %d barcode groups exported to database so far..." % i
                
        else:
            s.execute('update "%s" set "%s"="%s","%s"="%s","%s"="%s" where barcode="%s"' %
            (table, 'purity' + column, purities[i], 'fracClassified' + column, fracClassifiedLayer[i], 'mostCommon' + column, mostCommon[i], barcodes[i]))
            
            if (i % 1e3 == 0 and i != 0):
                sql.commit()        # Write reads to persistent DB
                print "Statistics for %d barcode groups exported to database so far..." % i
            
    sql.commit()        # Commit db changes to disk
    sql.close()         # Close the db
    
    print 'Phylogenetic statistics for column %s inserted into database!' % column

if __name__ == "__main__":
    
    ### Setup argument parser ###
     
    parser = argparse.ArgumentParser(description="""
    A script for analyzing phylogenetic data produced by Kraken.
    
    @author: Ben Demaree
    """, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('db', metavar='<db_filename>', type=str,
                        help='SQLite input database filename')
                        
    parser.add_argument('-t', metavar='table', type=str, default='kraken',
                        help='table in database to analyze [default: kraken]')
        
    args = parser.parse_args()      # Parse arguments
        
    sqliteDB = args.db      # DB filename
    table = args.t          # Table to analyze
        
    db = create_engine('sqlite:///' + sqliteDB)     # Create DB interface

    ### Compute phylogenetic statistics ###

    columns = ['domain', 'phylum', 'class', "order", 'family', 'genus', 'species']
    classifiedSizes = []    
    
    getNumClassified(table)
    
    for level in columns:
        # Compute phyologenetic statistics:
        purities, fracClassified, fracClassifiedLayer, mostCommon, barcodes = computeStats(table, level, classifiedSizes)
        
        # Insert data into database:
        insertStats(table, level, purities, fracClassified, fracClassifiedLayer, mostCommon, barcodes)
 
 
