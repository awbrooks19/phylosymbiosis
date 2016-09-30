#!/usr/bin/env python

### SCRIPT PARAMETERS ###
__author__ = "Andrew W. Brooks`"
__version__ = "1.0_utility_script"
__copyright__ = ""
__email__ = "andrew.w.brooks@vanderbilt.edu"
__status__ = "Template_Script"
### DATE: yy_mm_dd ###
__date__ = "16_7_13"

### IMPORT UTILITIES ###
import sys, getopt, os, csv, sets
from operator import itemgetter
from cogent.util.option_parsing import parse_command_line_parameters, make_option
from os import path
from os.path import join
import glob
import time

### IMPORT OPTIONAL ###
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt

#!!! RUNNING INFO !!!#
script_info = {}
script_info['brief_description'] = "Import Data Files from Microbiome Analysis"

#!!! DESCRIBE WHAT THE SCRIPT WILL DO HERE !!!#
script_info['script_description'] = ""

#!!! DESCRIBE SCRIPT OPTIONS FOR THE USER HERE !!!#
script_info['script_usage'] = [
("","Input a File:","%prog -i file.txt"),
("","Generate an Output Directory","%prog -o out_directory")]

script_info['output_description']= "X = [Headers, Data]"

#!!! ADD REQUIRED ARGUMENTS HERE !!!#
script_info['required_options'] = []

#!!! ADD OPTIONAL ARGUMENTS HERE !!!#
script_info['optional_options'] = [
    ### FOR INPUT FILE ###
    # make_option('-short_argument','--long_argument',default=None,type="existing_filepath",help='the input filepath'), 
    ### FOR OUTPUT DIRECTORY ###
    # make_option('-short_argument','--long_argument',default=None,type="new_dirpath",help='the output dirpath'), 
    
    make_option('-i','--input_file',default=None,type="existing_filepath",help='the biom table filepath'),
    make_option('-m','--map_file',default=None,type="existing_filepath",help='the map filepath'),
    make_option('-o','--output_file',default=None,type="new_dirpath",help='the output path for pdf'),
    ]
    
script_info['version'] = __version__

###################################################
### MAIN FUNCTION #################################
###################################################
#!! THIS IS WHERE THE MAIN PROGRAM OF THE SCRIPT RUNS
#!! YOU CAN PUT ANY CODE YOU WANT TO RUN IN THIS AREA
#!! MORE EXTENSIVE FUNCTIONS HOWEVER SHOULD BE CODED BELOW
def main():
    print "-._    _.--'\"''--._    _.--'\"''--._    _.--'\"''--._    _.--'\"''--._    _.--'\"''--._    _.--'\"''--._    _.--'\"''--._    _"
    print "-.  '-:'.'|'|\"':-.  '-:'.'|'|\"':-.  '-:'.'|'|\"':-.  '-:'.'|'|\"':-.  '-:'.'|'|\"':-.  '-:'.'|'|\"':-.  '-:'.'|'|\"':-.  '-:'"
    print " |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '."
    print " |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  "
    print " :_.' '.  '.:_ | :_.' '.  '.:_ | :_.' '.  '.:_ | :_.' '.  '.:_ | :_.' '.  '.:_ | :_.' '.  '.:_ | :_.' '.  '.:_ | :_.' '."
    print "-'       '-..,..-'       '-..,..-'       '-..,..-'       '-..,..-'       '-..,..-'       '-..,..-'       '-..,..-'      "
    print " 0===0                                                                                                            0===0 "
    print "  O=o                                                                                                              O=o  "
    print "   O                                                                                                                O   "
    print "  o=O                                             ANDREW W. BROOKS                                                 o=O  "
    print " 0===0                                     VANDERBILT GENETICS INSTITUTE                                          0===0 "
    print "  O=o                                    andrew.w.brooks(at)vanderbilt.edu                                         O=o  "
    print "   O                                            SETH R. BORDENSTEIN                                                 O   "
    print "  o=O                                                                                                              o=O  "
    print " 0===0                                                                                                            0===0 "
    print " -    -   -  - - --- - -  -   -    -"
    
    ### COGENT INITIALIZATION ###
    voption_parser, opts, args =\
        parse_command_line_parameters(**script_info)    
    #suppress_errors = opts.suppress_errors
    
    ### INPUT TABLE ARGUMENTS ###
    input_fps = []
    for input_fp in opts.input_file.split(','):
        input_fps.extend(glob.glob(input_fp))
        
    table_path = input_fps[0]
    print "TABLE FILEPATH: ", table_path
    
    ### INPUT MAP ARGUMENTS ###
    input_mps = []
    for map_fp in opts.map_file.split(','):
        input_mps.extend(glob.glob(map_fp))
        
    map_path = input_mps[0]
    print "MAP FILEPATH  : ", map_path
    
    ### INPUT OUTPUT ARGUMENTS ###
    output_fps = []
    for output_fp in opts.output_file.split(','):
        print output_fp
        output_fps.append(output_fp)
    print output_fps
    output_path = output_fps[0]
    print "OUTPUT FILEPATH: ", output_path
    
    ###################################################
    # Use BIOM Command Line Utilities to Add Map Data to Table File (tablepath.meta.biom)
    import os
    if not os.path.isfile(table_path + ".meta.biom"):
        print "biom add-metadata -i " + table_path + " -o " + table_path + ".meta.biom --sample-metadata-fp "+ map_path
        a = os.system("biom add-metadata -i " + table_path + " -o " + table_path + ".meta.biom --sample-metadata-fp "+ map_path)
        print a
    table_path = table_path + ".meta.biom"

    # Import Table
    from biom import load_table
    table = load_table(table_path) 
    
    get_table_info(table)
    observations = get_observations(table)
    samples = get_samples(table)
    metadata_categories = get_table_metadata_categories(table)
    metacats = get_table_metadata_categories(table)
    
    table_rel = get_table_relative(table)
    get_table_info(table_rel)
    #############
    ### INPUT ###
    phylum_idx = 2
    # Phylum_idx: 0 = Kingdom | 1 = Phylum | 2 = Class | 3 = Order | 4 = Family | 5 = Genus | 6 = Species
    #############
    
    table_rel_col = collapse_taxonomy(table_rel, phylum_idx)
    get_table_info(table_rel_col)
    
    table_rel_col_sort = table_rel_col.sort(axis='observation')
    get_table_info(table_rel_col_sort)
    
    taxa_in = table_rel_col_sort.ids(axis='observation')
    print 'Total Taxa: '+str(len(taxa_in))
    
    count = 0
    for values, id, metadata in table_rel_col_sort.iter(axis='sample'):
        print count, id, metadata['Both']
        count += 1
    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx

    # Figure Size
    fig = plt.figure(figsize=(120,60))
    ax = fig.add_subplot(121)

    # Get Range 
    taxa_range = np.arange(len(table_rel_col_sort.ids(axis='sample')))
    bottom = np.zeros(len(taxa_range))

    # Coloring
    cm = plt.get_cmap('rainbow') 
    cNorm  = colors.Normalize(vmin=0, vmax=len(table_rel_col_sort.ids(axis='observation')))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

    # Loop through Observations
    count = 0
    for values, id, metadata in table_rel_col_sort.iter(axis='observation'):
        if count % 25 is 0:  print str(count), " out of ", str(len(observations)), " Bacterial Taxon Plotted"
        values = np.array(values)
        colorval = scalarMap.to_rgba(count)
        plt.bar(taxa_range, values, color=colorval, bottom=bottom, label=samples, linewidth=0)
        bottom = np.add(bottom, values)
        count += 1

    axes = plt.gca()
    axes.set_xlim([0,len(taxa_range)])
    axes.set_ylim([0,1])
    plt.title('Taxonomy Bar Chart', fontsize=50)
    plt.ylabel('Class Relative Abundance', fontsize=40)
    plt.xlabel('Sample', fontsize=40)
    plt.legend(taxa_in, loc="upper right", bbox_to_anchor=(1.5,1))
    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()
    
    return


# Print info about table
def get_table_info(table_in):
    # Observation info
    print 'Total Observations: '+str(len(table_in.ids(axis='observation')))
    # Sample info
    print 'Total Samples: '+str(len(table_in.ids(axis='sample')))
    # Table info
    print 'Total Counts: '+str(table_in.sum())
    print 'Non-Zero Entries: '+str(table_in.nnz)
    print 'Table Density: '+str(table_in.get_table_density())

# Return List of Observations
def get_observations(table_in):
    obs = table_in.ids(axis='observation')
    return obs

# Return List of Samples
def get_samples(table_in):
    samp = table_in.ids(axis='sample')
    return samp

# Get Metadata Categories
def get_table_metadata_categories(table_in):
    [(values, id, metadata) for values, id, metadata in table_in.iter()]
    metac = []
    for md in metadata.keys():
        metac.append(str(md[:]))
    return metac

# Return Table as Relative Abundance
def get_table_relative(table_in):
    tablerel = table_in.norm(axis='sample', inplace=False)
    return tablerel 
    
# Return numpy 2D array of table counts
import numpy as np
def get_table_counts(table_in):
    # Get table shape
    sh = table_in.shape
    # Create Empty Numpy Table Filled with Zeros
    ztable = np.zeros((sh[1],sh[0]), dtype=np.int)
    # Add each Samples Counts
    cur_sample = 0
    cur_sum = 0
    for values, id, metadata in table.iter(axis='sample'):
        ztable[cur_sample] = values
        cur_sample += 1
    return ztable

# Collapse Taxonomy Function
def collapse_taxonomy(table, tax_level):
    # Collapse Function
    collapse_f = lambda id_, md: '; '.join(md['taxonomy'][:tax_level + 1])
    # Generate & Return Collapsed Table
    tabletaxacollapse = table.collapse(collapse_f, axis='observation',norm=False)
    return tabletaxacollapse


###################################################
### CRITICAL PIPELINES ############################
###################################################

### PIPELINE TO INPUT FILES ###
# This function should serve as a foundation for importing a wide variety of filetypes.
# The filepath should be provided and an integer indicating the type of file.
# Basic input should be tested for that file type, and the quality of the input
# should be validated.
def pipe_in(fpIn, typeIn=0):
    
    return
    
### PIPELINE FOR PANDAS INPUT ### 
# This function / class should carry out the input of pandas data frames and 
# then function as a user interface to display and clean data.
class pandas_class():
    
    ##########################
    ### PANDAS Constructor ###
    def __init__(self):
        self.df = None
        return
    
    ### PANDAS MAKE EMPTY DF ### construct a dataframe given lists of colsIn X indexIn
    def empty(self, colsIn, indexIn):
        self.df = pd.DataFrame(columns=colsIn, index=indexIn)
        
    ### PANDAS DISPLAY OPTIONS ### set display options for pandas dataframes
    def display_options(self, textWrap = False, maxRows = 10, maxCols = 10):
        pd.set_option('expand_frame_repr', textWrap)   # SET TEXT WRAPPING
        # pd.get_option("display.max_rows")            # GET MAXIMUM ROWS DISPLAYED
        pd.set_option("display.max_rows", maxRows)     # SET MAXIMUM ROWS DISPLAYED
        # pd.get_option("display.max_columns")         # GET MAXIMUM COLUMNS DISPLAYED
        pd.set_option("display.max_columns", maxCols)  # SET MAXIMUM COLUMNS DISPLAYED
        return
    
    ### PANDAS PIPELINE ### user interface for inputting pandas dataframe
    def input(fpIn, sepIn='\t', indexCol=None, skipRows=0):
        ### SET DEFAULT OPTIONS
        nRowIn = 10       # SET NUMBER OF ROWS TO DISPLAY
        nColIn = 10       # SET NUMBER OF COLUMNS TO DISPLAY
        verboseIn=False   # SET PANDAS VERBOSE MODE
        inputInt = 1000   # USER INTERFACE VARIABLE (INT)
        ### ENTER LOOP FOR USER INTERFACE...
        while inputInt != 0:
            #######################
            ### OPTIONS ###########
            ### 1 - INDEX COLUMN
            if inputInt == 1: boolInt, indexCol = input_int(qStr = "      --- SELECT INDEX COLUMN ---> ")
            ### 2 - SKIP ROWS
            if inputInt == 2: boolInt, skipRows = input_int(qStr = "      --- SKIP # ROWS ---> ")
        
            ### 3 - 
            #######################
            ### INPUT DATAFRAME
            dfIn = pd.read_csv(fpIn, sep=sepIn, index_col=indexCol, skiprows=skipRows, verbose=verboseIn)
            pd_display_options(maxRows = nRowIn, maxCols = nColIn, textWrap = False)
            ### DISPLAY TABLE
            print dfIn
            #######################
            ### PRINT USER OPTIONS ###
            print "   0. CONTINUE - Accept Current Table"
            print "   1. INDEX COLUMN - Select Index Column"
            print "   2. SKIP ROWS - Skip First X Rows"
            boolInt, inputInt = input_int(qStr = "   --- INPUT OPTION ---> ")
    
        return dfIn
    
    ### PANDAS FILL ROW BY INDEX ### fills data into indexIn row
    # dataIn must fill all columns
    def row_fill(self, iIn, dataIn):
        self.df.loc[iIn] = dataIn
    
    ### PANDAS FILL COLUMN  ### fills data into colIn column
    # dataIn must fill all index rows
    def col_fill(self, cIn, dataIn):
        self.df[cIn] = dataIn
    
    ### PANDAS FILL BOTH COLUMN AND INDEX POSITION  ### fills data specific location
    # dataIn must be a single point (i.e. integer or string)
    def pos_fill(self, iIn, cIn, dataIn):
        self.df.loc[iIn, cIn] = dataIn
    
    ### PANDAS EXPORT TO CSV  ### exports df as csv to locIn
    def export_csv(self, locIn):
        self.df.to_csv(locIn)

###################################################
### ADDITIONAL FUNCTIONS ##########################
###################################################

###################################################
### INPUT #########################################

### INPUT FUNCTION ### import delimited file
def read_csv(fpIn):
    csvIn = csv.reader(open(fpIn), delimiter='\t') # OPEN FILEPATH
    count = 0                                      # TRACK LINE NUMBER
    for rowIn in csvIn:                            # LOOP THROUGH ROWS
        # OPTIONAL WAY TO HANDLE ###
        if count <= 1:
            t1 = pd.Series(rowIn)
            print t1
        else: break
        #############################
        count += 1                                 # INCREMENT LINE NUMBER 
    return

###################################################
### OUTPUT #########################################

### OUTPUT FUNCTION ### enumerate through list and print #. value
def print_enumlist(listIn):
    for eIdx, eVal in enumerate(listIn): print "  "+str(eIdx)+". "+str(eVal)

###################################################
### DATE & TIME ###################################

### DATE & TIME FUNCTION ### return date & time as user readable string (returns string)
def date_print():
    dateStr = "Date: "+date_date()+" Time: "+date_time()
    return dateStr

### DATE & TIME FUNCTION ### get date (returns month_day_year string) 
def date_date():
    return time.strftime("%m_%d_%Y")

### DATE & TIME FUNCTION ### get time (return hour_minute_second string)
def date_time():
    return time.strftime("%H:%M:%S")

#################################################
### DIRECTORY ###################################

### DIRECTORY PIPELINE ### make output directory (returns bool and path string)
def dir_pipe(dir_path, overwrite=False):
    dir_path = dir_backslash(dir_path)           # Append backslash to end
    if not dir_check(dir_path):                  # If directory does not exist...
        dir_make(dir_path)                       # Then make directory
    elif overwrite:                              # Else if does exist and overwrite...
        dir_remove(dir_path);dir_make(dir_path)  # Then remove and remake
    if dir_check(dir_path): return True, dir_path# If Exists now then return
    else: return False, ""                       # Else return false

### DIRECTORY FUNCTION ### check if directory exists (returns bool true==exists)
def dir_check(dirPath):
    if os.path.isdir(dirPath): return True

### DIRECTORY FUNCTION ### make_directory (returns string of cwd/dir_path)
def dir_make(dirPath):
    os.makedirs(dirPath)
    
### DIRECTORY FUNCTION ### remove_directory (returns string of cwd/dir_path)
def dir_remove(dirPath):
    shutil.rmtree(dirPath)

### DIRECTORY FUNCTION ### append '/' to dir_path if not present (return string dir_path)
def dir_backslash(dirPath):
    if dirPath[-1] != '/': dirPath = dirPath+'/' # ADD '/' TO END IF NOT PRESENT
    return dirPath

### DIRECTORY FUNCTION ### get fullpath to directory (returns string of cwd/dir_path)
def dir_fullpath(dirPath):
    return os.path.abspath(dirPath)
    
### DIRECTORY FUNCTION ### get fullpath to cwd directory (returns string of cwd/dir_path)
def dir_localpath(dirPath):
    return str(os.getcwd())+'/'+dirPath # ADD FULLPATH

#################################################
### FILE ########################################

### FILE FUNCTION ### check if file exists (return bool)
def file_check(filePath):
    return os.path.isfile(filePath)

### FILE FUNCTION ### get seconds since last access of file (number)
def file_lasttime(filePath):
    return os.path.getatime(filePath)

### FILE FUNCTION ### get after last '/' in filepath as string
def file_basename(filePath):
    return os.path.basename(filePath)
    
#################################################
### USER INPUT ##################################

### USER INPUT FUNCTION ### get raw user input, evaluate if int (return bool, int)
def user_input_int(qStr = "  --- INPUT INTEGER ---> "):
    inputIn = raw_input(qStr)
    return eval_int(inputIn)
    
### USER INPUT FUNCTION ### get raw user input, evaluate if int range (return bool, [start,finish])
def user_input_ints(minRange=0,maxRange=10000000):
    inputIn = raw_input("  --- Input Integer Range(s): ex. 1:3,4,7: (: for all) ---> "); indexOut = []
    inputIn = inputIn.split(',')           # Split by ','
    for inSx, inS in enumerate(inputIn):   # for each split...
        inS2 = inS.split(':'); sectStore=[]# Split by ':'
        for inS3x,inS3 in enumerate(inS2): # for each split...
            eBool, eVal = eval_int(inS3)   # Evaluate if int
            if not eBool:                  # if not int...
                if inS3 == '':             # if ' ' or ':' then...
                    if inS3x == 0: eVal = minRange; eBool=True # if first then minRange
                    else: eVal = maxRange; eBool=True  # if last then maxRange
                else: print "  Not an Integer (Ignored): "+str(inS3); break
            sectStore.append(eVal)
        if eBool == True: indexOut.append(sectStore)
    return indexOut

#################################################
### EVALUATE ####################################

### EVALUATE FUNCTION ### returns true, # if number can be coerced to int
def eval_int(intIn):
    try: int(intIn); return True, int(intIn)     # Try coercing to int - return true, int
    except ValueError: return False, None        # Otherwise return false, None

    
if __name__ == "__main__":
    main()

