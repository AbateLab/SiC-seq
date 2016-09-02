# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 13:40:55 2016

@author: Ben
"""

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

if __name__ == "__main__":
    
    tax_levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    
    ### Setup argument parser ###
     
    parser = argparse.ArgumentParser(description="""
    A script for converting a kraken output file and generating a krona chart.
    
    The script generates two charts: one for all reads pooled together and
    another for barcode groups, classified by the specified taxonomic level.
    
    Assumes Kraken output file is in MPA format and krakenAnalysis.py has been
    run on kraken table in db.

    ***Dependencies***
    kronaTools, sqlite3
    
    @author: Ben Demaree
    """, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('db', metavar='<db_filename>', type=str,
                        help='SQLite input database filename')
    
    parser.add_argument('-k', metavar='kraken_out', type=str, default='out.kraken',
                        help='Kraken output filename [default: out.kraken]')
    
    parser.add_argument('-t', metavar='kraken_table', type=str, default='kraken',
                        help='database table containing Kraken data [default: kraken]')
                            
    parser.add_argument('-l', metavar='tax_level', type=str, 
                        choices = tax_levels, default='genus',
                        help='lowest taxonomic level to include in chart [default: genus]')
                        
    parser.add_argument('-f', metavar='frac_classified', type=float, default=0.0,
                        help='minimum fraction of classified reads in a bargroup to include in chart [default: 0]')
                        
    parser.add_argument('-n', metavar='center_label', type=str, default='all',
                        help='label for center of chart [default: all]')
        
    args = parser.parse_args()      # Parse arguments
    
    sqliteDB = args.db              # DB filename
    tableK = args.t                 # Table for Kraken data
    krakenFile = args.k             # Kraken output filename
    fracClass = args.f              # Krona output chart name
    level = args.l                  # Krona output chart name
    centerLabel = args.n            # Label for center of Krona chart
    
    db = create_engine('sqlite:///' + sqliteDB)     # Create DB interface

    kronaChart = 'kronaChart_allreads.html'
    kronaChartBar = 'kronaChartBar_by%s_min%0.1f%%classfied.html' % (level, fracClass * 100)             
    
    outFileName = 'kronaOut.txt'
    outFileNameBar = 'kronaOutBar.txt'
    
    ### Generate charts ###
    
    # Chart for all reads:    
    
    krakenFile = open(krakenFile, 'r')
    outFile = open(outFileName, 'w')
    
    orgDict = {}                            # Dictionary for storing taxonomic classifications by level
    levelInd = tax_levels.index(level)      # Taxonomic level index

    for line in krakenFile:
        
        if line.split('\t')[1][:4] == 'root':       # Ignore reads classified as 'root'
            continue
        
        else:
            phylo = line.split('\t')[1][3:]         # Remove bargroup label
            phylo = re.sub('\|\w__', '\t', phylo)   # Remove phylogenetic delimiters
            phylo = phylo.replace('_', ' ')         # Replace underscores with spaces
            phylo = phylo.replace('\n', '')         # Remove newline
            
            outFile.write(phylo + '\n')        # Write to Krona chart file

            try:
                if phylo.split('\t')[levelInd] not in orgDict:
                    phyloTruncated = phylo.split('\t')[:levelInd + 1]
                    orgDict[phylo.split('\t')[levelInd]] = '\t'.join(phyloTruncated)
                    
            except IndexError:
                continue
     
    outFile.close()
    krakenFile.close() 
    
    makeChart = subprocess.call('perl /home/ben/programs/Krona/KronaTools/scripts/ImportText.pl -q %s -o %s -n %s' % (outFileName, kronaChart, centerLabel), shell = True)

    # Chart by barcode:
    
    df = pd.read_sql_query('select "%s" from "%s" where fracClassified>"%s" and "%s" is not null group by barcode' % 
                            ('mostCommon' + level, tableK, fracClass, 'mostCommon' + level), db)    
    
    mostCommon = df['mostCommon' + level]       # List of most common entry by taxonomic level, grouped by barcode
    
    outFileBar = open(outFileNameBar, 'w') 
    
    for org in mostCommon:
        outFileBar.write(orgDict[org.replace('_', ' ')] + '\n')          # Write taxonomic classification of bargroup to file
        
    outFileBar.close()
    
    makeChart = subprocess.call('perl /home/ben/programs/Krona/KronaTools/scripts/ImportText.pl -q %s -o %s -n %s' % (outFileNameBar, kronaChartBar, centerLabel), shell = True)
    
