# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 14:16:05 2016

@author: Ben Demaree
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
import re

### Pandas imports
import pandas as pd
from sqlalchemy import create_engine

# To list available fonts:
#import matplotlib.font_manager
#list = matplotlib.font_manager.get_fontconfig_fonts()
#names = [matplotlib.font_manager.FontProperties(fname=fname).get_name() for fname in list]
#print names

############## Functions ##############

Qscore = dict((chr(i),i-33) for i in range(33,90))      # Create quality score dictionary
    
def writeFASTA (df, outFile):
# Writes read ID and sequence info to a FASTA file for blast analysis

    FASTAfile = open(outFile, "w")
    
    for i in range(len(df['barcode'])):
        
        ID = str(df['barcode'][i]) + '-' + str(df['barID'][i])
        
        FASTAfile.write(">" + ID + '-1\n')
        FASTAfile.write(df['seq1'][i] + '\n')
        
        if not SE:
            FASTAfile.write('>' + ID + '-2\n')
            FASTAfile.write(df['seq2'][i] + '\n')
        
    FASTAfile.close()
    
def krakenTable (table, tableK):
# Creates a new table in the database for storing kraken phylogenetic data
    
    sql = sqlite3.connect(sqliteDB, isolation_level = 'Exclusive')
    s = sql.cursor()    # Create cursor
    
    tabExists = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table' AND name='%s'" % tableK, db)
    
    if len(tabExists['name']) > 0:
        overwrite = raw_input('%s table already exists in this database. Overwrite? (Y/N) ' % tableK)
        
        yes = set(['yes', 'y', 'Y'])
        no = set(['no','n', 'N'])
    
        if overwrite in yes:
            s.execute('drop table "%s"' % tableK)
        elif overwrite in no:
            print "Please change table names and run script again."
            raise SystemExit
        else:
            print "Please respond with 'Y' or 'N'."
            raise SystemExit
    
    try:
        print 'Creating table for Kraken data...'
            
        # Create table for kraken data:
        if not SE:    
            s.execute('create table "%s" as select barcode,barID,barSize,seq1,seq2,diversity from "%s"' % (tableK, table))
        else:
            s.execute('create table "%s" as select barcode,barID,barSize,seq1,diversity from "%s"' % (tableK, table))
    except sqlite3.Error:
        pass
    
    try:    
        # Add columns for phylogenetic data:
        s.execute('alter table "%s" add column classified TEXT default "N"' % tableK)  # Classification status (Y or N)
        s.execute('alter table "%s" add column domain TEXT' % tableK)          # Domain/superkingdom
        s.execute('alter table "%s" add column phylum TEXT' % tableK)          # Phylum
        s.execute('alter table "%s" add column class TEXT' % tableK)           # Class
        s.execute('alter table "%s" add column "order" TEXT' % tableK)         # Order
        s.execute('alter table "%s" add column family TEXT' % tableK)          # Family
        s.execute('alter table "%s" add column genus TEXT' % tableK)           # Genus
        s.execute('alter table "%s" add column species TEXT' % tableK)         # Species   
    except sqlite3.Error:
        pass
    
    try:       
        # Create index on barcodes:
        s.execute('create index "%s" on "%s" (barcode)' % (tableK + 'barIndex', tableK))
    except sqlite3.Error:
        pass

    try:       
        # Create index on barID:
        s.execute('create index "%s" on "%s" (barID)' % (tableK + 'barIndex', tableK))
    except sqlite3.Error:
        pass
    
    print 'Table %s set up successfully.' % tableK    
    
    sql.commit()        # Commit db changes to disk
    sql.close()         # Close the db
    
def kraken (inFile = 'kraken.fasta', tempFile = 'temp.kraken', outFile = 'out.kraken'):
# Runs kraken and parses the output data for export into the sqlite database    
    
    # Run kraken and translate results into readable file:
    krakenRun = """kraken --db %s %s --threads 6 --quick --min-hits 2 --output %s""" % (krakenDB, inFile, tempFile)
    krakenTrans = """kraken-translate --db %s %s --mpa-format > %s""" % (krakenDB, tempFile, outFile)
    
    krakenExec = subprocess.call(krakenRun, shell = True)
    krakenExec = subprocess.call(krakenTrans, shell = True)    

    # Parse output file from kraken-translate:

    sql = sqlite3.connect(sqliteDB, isolation_level = 'Exclusive')
    s = sql.cursor()    # Create cursor
    # These settings allow DB to be stored in memory and written periodically to the disk: 
    s.execute('PRAGMA synchronous = 0')
    s.execute('PRAGMA journal_mode = OFF')
            
    # Insert purity and most common species information into database:
    
    out = open(outFile, 'r')        # Open kraken output for reading
    
    count = 0           # Number of reads examined
    rootCount = 0       # Number of reads with 'root' classification (discarded)
    
    barcodes_ID = []    # Reads already analyzed (skip second mate in a read-pair)
    
    for line in out:
        
        ID = line.split('\t')[0].split('-')
        org = line.split('\t')[1]
        
        if org[:4] == 'root':    # Ignore reads with no phylogenetic classifiers
            rootCount += 1
            continue
        
        barcode = ID[0]     # Extract barcode
        barID = ID[1]       # Extract barcode ID
        
        # Skip second mate if a read pair if already classified
        if not SE:
            
            if barcode + barID in barcodes_ID:
                continue
            
            else:
                barcodes_ID += [barcode + barID]
       
        count += 1  # Increment read counter
        
        if (count % 1e4 == 0 and count != 0):
            sql.commit()        # Write reads to persistent DB
            print "Phylogenetic data from %d reads inserted into database so far..." % count
       
        # Clean up formatting:
        org = re.sub('\w__', '', org).replace('\n', '').split('|')
        
        # Parse taxonomy into levels:
                    
        # Set default classifications to empty string:
        phylum = ''; Class = ''; order = ''; family = ''; genus = ''; species = ''
        
        domain = org[0]
        
        try:        
            phylum = org[1]
        except IndexError:
            s.execute('update "%s" set classified="Y",domain="%s",phylum="%s",class="%s","order"="%s",family="%s",genus="%s",species="%s" where barcode="%s" and barID="%s"' %
            (tableK, domain, phylum, Class, order, family, genus, species, barcode, barID))
            continue

        try:        
            Class = org[2]
        except IndexError:
            s.execute('update "%s" set classified="Y",domain="%s",phylum="%s",class="%s","order"="%s",family="%s",genus="%s",species="%s" where barcode="%s" and barID="%s"' %
            (tableK, domain, phylum, Class, order, family, genus, species, barcode, barID))
            continue
            
        try:
            order = org[3]
        except IndexError:
            s.execute('update "%s" set classified="Y",domain="%s",phylum="%s",class="%s","order"="%s",family="%s",genus="%s",species="%s" where barcode="%s" and barID="%s"' %
            (tableK, domain, phylum, Class, order, family, genus, species, barcode, barID))
            continue
        
        try:        
            family = org[4]
        except IndexError:
            s.execute('update "%s" set classified="Y",domain="%s",phylum="%s",class="%s","order"="%s",family="%s",genus="%s",species="%s" where barcode="%s" and barID="%s"' %
            (tableK, domain, phylum, Class, order, family, genus, species, barcode, barID))
            continue
        
        try:
            genus = org[5]
        except IndexError:
            s.execute('update "%s" set classified="Y",domain="%s",phylum="%s",class="%s","order"="%s",family="%s",genus="%s",species="%s" where barcode="%s" and barID="%s"' %
            (tableK, domain, phylum, Class, order, family, genus, species, barcode, barID))
            continue
        
        try:
            species = org[6]
            s.execute('update "%s" set classified="Y",domain="%s",phylum="%s",class="%s","order"="%s",family="%s",genus="%s",species="%s" where barcode="%s" and barID="%s"' %
            (tableK, domain, phylum, Class, order, family, genus, species, barcode, barID))
            continue
        except IndexError:
            s.execute('update "%s" set classified="Y",domain="%s",phylum="%s",class="%s","order"="%s",family="%s",genus="%s",species="%s" where barcode="%s" and barID="%s"' %
            (tableK, domain, phylum, Class, order, family, genus, species, barcode, barID))
            continue

    sql.commit()        # Commit db changes to disk one last time
    sql.close()         # Close the db
    
    print 'Phylogenetic data from %d reads inserted into database!' % count
    print '%d reads had no classifiable LCA and were discarded.' % rootCount    
    
    out.close()

if __name__ == "__main__":
    
    ### Setup argument parser ###
     
    parser = argparse.ArgumentParser(description="""
    A script for running the Kraken phylogenetic classifier and importing data 
    into an SQLite database.
    
    You must manually set the desired database to use in the script.

    ***Dependencies***
    sqlite3, kraken
    
    @author: Ben Demaree
    """, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('db', metavar='<db_filename>', type=str,
                        help='SQLite input database filename')
                        
    parser.add_argument('-t', metavar='reads_table', type=str, default='reads',
                        help='table in database to analyze [default: reads]')
                        
    parser.add_argument('-k', metavar='kraken_table', type=str, default='kraken',
                        help='new table name containing Kraken data [default: kraken]')
    
    parser.add_argument('--F', action='store_true', default=False,
                        help='option to use existing kraken.fasta file')
                                                                     
    parser.add_argument('--SE', action='store_true', default=False,
                        help='flag for a single-end read analysis') 
        
    args = parser.parse_args()      # Parse arguments
        
    sqliteDB = args.db      # DB filename
    table = args.t          # Table to analyze
    tableK = args.k         # Table for Kraken data
    F = args.F              # Option to use existing kraken.fasta file
    SE = args.SE            # Option to analyze a single-end read
        
    db = create_engine('sqlite:///' + sqliteDB)     # Create DB interface
    
    ### Kraken analysis ###
    
    krakenDB = '/drive1/genomes/KrakenDBs/minikraken_20141208'
    
    krakenTable(table, tableK)  # Create Kraken table
    
    if not SE:
        df = pd.read_sql_query('select barcode,barID,seq1,seq2 from "%s"' % table, db)
    else:
        df = pd.read_sql_query('select barcode,barID,seq1 from "%s"' % table, db)
   
   # Write reads to kraken.fasta (if selected):
    if not F:
        writeFASTA(df, 'kraken.fasta')
        
    kraken()    # Run Kraken and insert data into sqlite db
    
