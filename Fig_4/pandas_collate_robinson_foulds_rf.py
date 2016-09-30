#!/usr/bin/env python
# Made on 8/27/15

### SCRIPT PARAMETERS
__author__ = "Andrew Brooks`"
__version__ = "1.0_brooks_utils"
__copyright__ = ""
__email__ = "andrew.w.brooks@vanderbilt.edu"
__status__ = ""

### IMPORT UTILITIES ###
import sys, getopt, os, csv, sets, itertools                                      # General Python Utility Packages
from cogent.util.option_parsing import parse_command_line_parameters, make_option # COGENT script automation
import glob                                                                       # To perform REGEX searches
import time                                                                       # To get date and time
import shutil                                                                     # To Remove Directory

#################################################
############## COGENT RUNNING INFO ##############
#################################################

script_info = {}
script_info['brief_description'] = "Bioinformatics Utility Script by Andrew Brooks"

# DISPLAYED SCRIPT DESCRIPTION
script_info['script_description'] = \
"//////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n"\
"####################################### SCRIPT INFO ########################################\n"\
"                       Bioinformatics Utility Script by Andrew Brooks                       \n"\
"                         - awbrooks19(at)gmail.com                                          \n"\
"                         - Vanderbilt: Center for Human Genetics                            \n"\
"                                                                                            \n"\
"           ----------------            DESCRIPTION               ----------------           \n"\
"                                                                                            \n"\
"           ----------------             WARNINGS                 ----------------           \n"\
"  - Scripts are not written to catch all possible user errors.                              \n"\
"    Use wisely, consiously, and at your own discretion!                                     \n"\
"                                                                                            \n"\
"  - The [-f] option will overwrite the current output folder, use at your own risk!         \n"\
"                                                                                            \n"

# GUIDE HOW TO USE SCRIPT OPTIONS
script_info['script_usage'] = [\
 ("Run With User Input","Add -v","%prog -v"),
 ("Import from single file","provide path to single file","%prog -i in.txt"),
 ("Import from list of files","provide comma-separated list of files","%prog -i in1.txt,in2.txt"),
 ("Import from all files with regex","provide regex expression to filter files","%prog -i \"*.txt\"")]

script_info['output_description']= "OPTIONAL"

### REQUIRED OPTIONS ###
script_info['required_options'] = [
 # INPUT [-i] FILEPATH
 # make_option('-i','--input_fps',help='the input filepaths'),
]

### OPTIONAL OPTIONS ###
script_info['optional_options'] = [

 # OUTPUT [-o] FILEPATH #
 make_option('-o','--output_fp',help='Output Directory'),
 
 # INPUT [-i] FILEPATH(S) #
 make_option('-i','--input_fps',help='Input Filepaths'),
 
 # OVERWRITE OUTPUT [-o] #
 make_option('-f','--overwrite',action='store_true', help='Overwite the Output Directory [default: %default]',default=False),
]

script_info['version'] = __version__

#################################################
################ MAIN FUNCTION ##################
#################################################
def main():

    # INPUT COMMAND LINE OPTIONS
    voption_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    ########### VOCAL ###########
    if opts.verbose:
        print script_info['script_description']
        print "####################################### INITIALIZATION #####################################"                                
        print date_print()
        print "Verbose: "+str(opts.verbose)
        print "PWD:    "+os.getcwd()
        
    ########## OUTPUT ###########
    # If output folder specified: make output folder
    if opts.output_fp:
        out_bool, out_path = dir_pipe(opts.output_fp, overwrite=opts.overwrite)
        if (out_bool == True) and (opts.verbose == True): print "Output: "+dir_fullpath(out_path) # print output
        elif opts.verbose == True: print "Output: Could not be created: "+out_path # if not found then 
    
    ### PRINT MODULES ###
    #ref_modules(justModuleIn=True, justFunctionIn=True)
    ### PRINT FUNCTIONS ###
    #ref_modules(justModuleIn=False, justFunctionIn=True)
    ### PRINT FUNCTIONS & DESCRIPTIONS ###
    #ref_modules(justModuleIn=False, justFunctionIn=False)
    
    ########## INPUT ############
    if opts.verbose: print "########################################## INPUT [-i] ######################################"
    input_input_fps = input_pipe(opts.input_fps, optTitle="[-i --input_fps]", verboseIn=opts.verbose)
    
    ######### CUSTOM ############
    if opts.verbose: print "########################################### CUSTOM #########################################"
    if input_input_fps: print_enumlist(input_input_fps)
    
    # GET PANDAS FILES IMPORTED # 
    pdInputs = []
    for inF in input_input_fps:
        pdInputs.append(pd_ui(inF))
    
    
    # STOCHASTIC FIRST, THEN HOST-MICROBE TREE COMPARISON #
    pd_summarize(pdInputs[0], locIn="Stochatic Comparison")
    pd_summarize(pdInputs[1], locIn="Host-Microbe Comparison")
    
    metricIn="R-F"

    hostMicrobeScore = pdInputs[1][metricIn][0]

    print 
    print "Host-Microbe Score:   " + str(hostMicrobeScore)
    print
    print "Better Score:         " + str(len(pdInputs[0][pdInputs[0][metricIn] < hostMicrobeScore]))
    print "Worse Score:          " + str(len(pdInputs[0][pdInputs[0][metricIn] > hostMicrobeScore]))
    print "Equiv Score:          " + str(len(pdInputs[0][pdInputs[0][metricIn] == hostMicrobeScore]))
    print "P-value better:       " + str(float(len(pdInputs[0][pdInputs[0][metricIn] < hostMicrobeScore]))/100000.0)
    print 
    print "Better\Equal Score:   " + str(len(pdInputs[0][pdInputs[0][metricIn] <= hostMicrobeScore]))
    print "Worse Score:          " + str(len(pdInputs[0][pdInputs[0][metricIn] > hostMicrobeScore]))
    print "P-value Better/Equal: " + str(float(len(pdInputs[0][pdInputs[0][metricIn] <= hostMicrobeScore]))/100000.0)
    print 
    print "Max Stochastic Metric:  " + str(max(pdInputs[0][metricIn]))
    
    #############################
    if opts.verbose:
        print "############################################ END ###########################################"
        print "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////////////"

##################################################################################################
##################################################################################################
############### CUSTOM FUNCTIONS #################################################################
##################################################################################################
##################################################################################################

"""      ### THESE ARE PACKAGES USED MOST OFTEN ###
### OPTIONAL UTILITIES ###
import numpy as np                    # NUMPY data storage and analysis
import scipy.stats as sp              # SCIPY statistical analyses
import statsmodels.api as sm          # STATSMODELS

### PLOTTING UTILITIES ###
import matplotlib.pyplot as plt       # MATPLOTLIB plotting tools
import seaborn as sns                 # SEABORN plotting tools
sns.set(color_codes=True)
"""      ##########################################
import pandas as pd                   # PANDAS dataframe tools

#################################################
################## PANDAS #######################!!!
### PANDAS PIPELINE ### import pandas dataframe and manipulate
def pd_ui(fpIn, sepIn='\t', indexCol=None, skipRows=0,
          verboseIn=False, textWrap=False):
    userInt = None
    pd.set_option('expand_frame_repr', textWrap)       # No Wrapping in display

    while userInt != 0:                     # While not option 0 (accept table)...
        dfIn = pd.read_csv(fpIn, sep=sepIn, index_col=indexCol, skiprows=skipRows,
                           verbose=verboseIn)
        print dfIn
        print "  ---------------------> "
        print "  0. Accept Table"           # Print Options
        print "  1. Set Delimiter    ['\\t']"
        print "  2. Set Index Column [None]"
        print "  3. Set Skip Rows    [0]   "
        userBool, userInt = input_int()
        if userBool:
            print userInt
            if userInt == 0: break
    return dfIn     # If option 0 (accept table) then return dataframe

def pd_summarize(dfIn, locIn=""):
    print "Summary: "+locIn
    print "  Number Columns: "+str(len(pd_columns(dfIn)))
    print "  Number Indices: "+str(len(pd_indices(dfIn)))
    return
"""
def pd_ui():
	# LABELS: String or [Strings,...]
	# AXIS: 0=Row | 1=Column
	# INPLACE: (True) acts on current dataframe and returns None
	# LEVEL: for nested dataframe
	dfIn.drop(labels, axis=0, inplace=True, level=None, errors='raise')

### PANDAS FUNCTION ### remove duplicate rows 
def pd_drop_duplicates():
    # KEEP: 'first', 'last', False, default 'first'
    # SUBSET: can provide list of only Row ID's to consider
	dfIn.drop_duplicates(subset, inplace=True, keep='first')
"""
def pd_values(dfIn):
    return dfIn.values
    
def pd_indices(dfIn):
    return dfIn.index
    
def pd_columns(dfIn):
    return dfIn.columns

def pd_describe(dfIn):
    print df.describe()
    
def pd_transpose(dfIn):
    return dfIn.T
    
def pd_sortindex_col(dfIn, ascendIn=True):
    return dfIn.sort_index(axis=1, ascending=ascendIn)

def pd_sortindex_row(dfIn, ascendIn=True):
    return dfIn.sort_index(axis=0, ascending=ascendIn)

def pd_daterange(startDate='20160101', numDays=5):
    return pd.date_range(startDate, periods=numDays)

##################################################################################################
##################################################################################################
############# GENERAL FUNCTIONS (Required) #######################################################
##################################################################################################
##################################################################################################

#################################################
################ REFERENCE ######################
### REFERENCE REFERENCES ### use tuple to organize functions
def ref_ref():
    optsRef=[]
    # ref_print()
    optsRef.append(ref_func(optNum = 0, optFunc = "ref_print", optHelp = "Print Out Function References for Module",
             optVars = [("optsRefIn","List(referenceTuples)","Input List of Function Reference Tuples"),
                        ("optsRefTitle=","String[None]","String of Module Title")],
             optRet = "None"))
    return optsRef

### REFERENCES FUNCTION ### print all reference modules
def ref_modules(justModuleIn=False, justFunctionIn=False):
    ref_print(ref_ref(), optsRefTitle="REFERENCES (ref_)", justModule=justModuleIn, justFunction=justFunctionIn)
    ref_print(ref_input(), optsRefTitle="INPUT (input_)", justModule=justModuleIn, justFunction=justFunctionIn)
    ref_print(ref_search(), optsRefTitle="SEARCH (search_)", justModule=justModuleIn, justFunction=justFunctionIn)

### REFERENCE FUNCTION ### print out references for type of functions
# optRef[0](numID) | optRef[1](functionName) | optRef[2](functionHelp) | optRef[3][list of input options] | optRef[4](outputType)
def ref_print(optsRefIn, optsRefTitle=None, justModule=False, justFunction=False):
    if optsRefTitle: print "------- "+optsRefTitle+" --------"
    if justModule == True: return
    for optRef in optsRefIn:
        print "  -"+str(optRef[0])+"- "+optRef[1]+"() : "+optRef[2]+" ---"
        if not justFunction:
            for funcInput in optRef[3]:
                if "=" in funcInput[0]: print "      > Optional: "+str(funcInput)
                else: print "      > Required: "+str(funcInput)
            print "      < Return: "+optRef[4]

### REFERENCES FUNCTION ### create ref option
def ref_func(optNum = 0, optFunc = "ex_func", optHelp = "ex help...",
             optVars = [("var1","type(var1)[Default]","var1 help..."),
                        ("var2=","type(var2)[Default]","var2 help...")],
             optRet = "return type..."):
    return (optNum,optFunc,optHelp,optVars,optRet)

#################################################
################## INPUT ########################
### REFERENCES INPUT ### use tuple to organize functions
def ref_input():
    optsRef=[]
    # input_pipe()
    optsRef.append(ref_func(optNum = 0, optFunc = "input_check", optHelp = "Split Input by ','",
             optVars = [("optIn","String(OptionInput)","Input Option as String to be Split"),
                        ("check_regex=","Boolean[False]","Perform REGEX Search on each Split"),
                        ("check_files=","Boolean[False]","Check if Files Exist at each Filepath")],
             optRet = "[ListOfValues]"))
    # input_check()
    optsRef.append(ref_func(optNum = 0, optFunc = "input_check", optHelp = "Split Input by ','",
             optVars = [("optIn","String(OptionInput)","Input Option as String to be Split"),
                        ("check_regex=","Boolean[False]","Perform REGEX Search on each Split"),
                        ("check_files=","Boolean[False]","Check if Files Exist at each Filepath")],
             optRet = "[ListOfValues]"))
    # input_pick()
    optsRef.append(ref_func(optNum = 1, optFunc = "input_pick", optHelp = "Pick Input Files to Use",
             optVars = [("inputIn","List[InputFilePaths]","List of Input Filepaths to Choose From")],
             optRet = "[ListOfChosenFilePaths]"))
    # input_int()
    optsRef.append(ref_func(optNum = 2, optFunc = "input_int", optHelp = "User Input of Int",
             optVars = [],
             optRet = "Boolean(isInt), int(intValue)"))
    ###
    return optsRef

##### FOR INPUT DATA THAT IS NOT FILENAMES #####
### INPUT STRING ###
# This function will take a Cogent opts... and parse it into list of strings... comma separated into a list of strings (or list of single string)
def input_string_data(optIn):
    optIn = optIn.split(',')
    return optIn

### INPUT FUNCTION ### get raw user input, evaluate if int (return bool, int)
def input_int_data(optIn):
    return eval_int(optIn)

##### FOR INPUT FILES AS PROGRAM ARGUMENTS ex. [-i] #####
### INPUT PIPELINE ### check if input option exists and find\select files
def input_pipe(optIn, optTitle="[-i --input_fps]", verboseIn=True):
    if optIn:
        inInputFps = input_check(optIn, check_files=True, check_regex=True)
        if verboseIn:                                     # if verbose...
            print "Input Files " + optTitle + ":  "+optIn # Print Input File regex
            print "Files Found:"
            print_enumlist(inInputFps)                  # Print all the input files...
            print "Select Files to Use:"
            inInputFps = input_pick(inInputFps)       # Pick Input to use
            print "Returned Files:"
            print_enumlist(inInputFps)
        return inInputFps
    else: print "Option not Found: "+optTitle; return False
    
### INPUT FUNCTION ### check if option must be split, search for all regex files, and if check_files then validate they exist
def input_check(optIn, check_regex=False, check_files=False):  # (returns list!)
    optIn = optIn.split(',')     # Split by ','
    if check_regex == True:      # If Regex...
        opt_regvals = []         # Storage regex finds
        for opt_val in optIn:    # For each input string...
            opt_regvals.extend(search(opt_val))   # Perform regex search and store
        optIn = opt_regvals      # Reset optIn to contain regex
    if check_files == True:      # If check...
        for opt_val in optIn:    # For each input string...
            if not os.path.exists(opt_val): print "File Missing: "+opt_val; optIn.remove(opt_val) # If not file or dir then print error
    return optIn   

### INPUT PIPELINE ### select input files to keep (user input, returns list)
def input_pick(inputIn):
    intsPick = input_ints(minRange=0,maxRange=len(inputIn)-1); inputOut=[]
    for intR in intsPick:
        if len(intR) == 1: inputOut.append(inputIn[intR[0]]) # validate in range and store value        
        else: 
            for cApp in inputIn[intR[0]:intR[1]+1]: inputOut.append(cApp)
    return inputOut

##### WITHIN PROGRAM QUERIES TO USER #####
### INPUT FUNCTION ### get raw user input, evaluate if int (return bool, int)
def input_int():
    inputIn = raw_input("  --- Input Integer ---> ")
    return eval_int(inputIn)

### INPUT FUNCTION ### get raw user input, evaluate if int range (return bool, [start,finish])
def input_ints(minRange=0,maxRange=10000000):
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

### INPUT FUNCTION ### get raw user input, returned as string, trailing white trimmed
def input_user(inQuery):
    return raw_input(inQuery)

#################################################
################# SEARCH ########################
### REFERENCE SEARCH ###
def ref_search():
    optsRef=[]
    # search()
    optsRef.append(ref_func(optNum = 0, optFunc = "search", optHelp = "Perform REGEX Search of Input String",
             optVars = [("searchStr","String(filepath)","Input String Filepath to REGEX Search")],
             optRet = "List[String(Results)]"))
    return optsRef
    
### SEARCH FUNCTION ### regex search files (returns list of files matching regex expression)
def search(searchStr):
    return glob.glob(searchStr)

#################################################
################ DIRECTORY ######################
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
################## FILE #########################
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
################## DATE #########################

### DATE FUNCTION ### return date & time as user readable string (returns string)
def date_print():
    dateStr = "Date: "+date_date()+" Time: "+date_time()
    return dateStr

### DATE FUNCTION ### get date (returns month_day_year string) 
def date_date():
	return time.strftime("%m_%d_%Y")

### DATE FUNCTION ### get time (return hour_minute_second string)
def date_time():
	return time.strftime("%H:%M:%S")

#################################################
################ EVALUATE #######################
### EVALUATE FUNCTION ### returns true, # if number can be coerced to int
def eval_int(intIn):
    try: int(intIn); return True, int(intIn)     # Try coercing to int - return true, int
    except ValueError: return False, None        # Otherwise return false, None

### EVALUATE FUNCTION ### converts number [0:255] into ASCII character
def eval_numascii(numIn):
    if numIn >= 0 and numIn <256: return chr(numIn)
    else: print "  Number out of ASCII Range [0:255]: "+str(numIn); return

#################################################
################## PRINT ########################
### PRINT FUNCTION ### enumerate through list and print #. value
def print_enumlist(listIn):
    for eIdx, eVal in enumerate(listIn): print "  "+str(eIdx)+". "+str(eVal)

#################################################
################## NUMBERS ######################
### NUMBERS FUNCTION ### get absolute number
def num_abs(numIn):
    return abs(numIn)

### NUMBERS FUNCTION ### takes in a list or array and returns generator (for i in generator:)
def num_chunks(listIn, chunkSize):
    for chunkBlock in range(0, len(listIn), chunkSize): 
        yield listIn[chunkBlock:chunkBlock+chunkSize]

#################################################
#################### END ########################
#################################################
       
if __name__ == "__main__":
    main()