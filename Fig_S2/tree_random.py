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

### OPTIONAL UTILITIES ###
import pandas as pd                   # PANDAS dataframe tools
import numpy as np                    # NUMPY data storage and analysis
import scipy.stats as sp              # SCIPY statistical analyses
import statsmodels.api as sm          # STATSMODELS

### PLOTTING UTILITIES ###
import matplotlib.pyplot as plt       # MATPLOTLIB plotting tools
import seaborn as sns                 # SEABORN plotting tools
sns.set(color_codes=True)


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
"  Generates a number of trees with random shuffled topologies for the given leaves.         \n"\
"  These trees were used as a null hypothesis of stochastic sample inter-relationships       \n"\
"  for the Robinson-Foulds Test.                                                             \n"\
"  The branch distances are '1' for all, and all splits ar bifurcations.                     \n"\
" \n"\
"           ----------------             WARNINGS                 ----------------           \n"\
"  - Scripts are not written to catch all possible user errors.                              \n"\
"    Use wisely, consiously, and at your own discretion!                                     \n"\
"                                                                                            \n"\
"  - The [-f] option will overwrite the current output folder, use at your own risk!         \n"\
"                                                                                            \n"

# GUIDE HOW TO USE SCRIPT OPTIONS
script_info['script_usage'] = [\
 #("Run With User Input","Add -v","%prog -v"),
 #("Import from single file","provide path to single file","%prog -i in.txt"),
 #("Import from list of files","provide comma-separated list of files","%prog -i in1.txt,in2.txt"),
 #("Import from all files with regex","provide regex expression to filter files","%prog -i \"*.txt\""),
 ]

script_info['output_description']= "OPTIONAL"

### REQUIRED OPTIONS ###
script_info['required_options'] = [
 # INPUT [-i] FILEPATH
 # make_option('-i','--input_fps',help='the input filepaths'),
 
 
 ##### CUSTOM #####
 make_option('-n','--num_trees',help='the number of trees to generate'),
 make_option('-i','--leaf_names',help='the names (comma separated) of leaves ex. a,b,c'),
 
]

### OPTIONAL OPTIONS ###
script_info['optional_options'] = [

 # OUTPUT [-o] FILEPATH #
 make_option('-o','--output_fp',help='Output Directory'),
 
 
 # INPUT [-i] FILEPATH(S) #
 #make_option('-i','--input_fps',help='the input filepaths'),
 
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
    #if opts.verbose: print "########################################## INPUT [-i] ######################################"
    #input_input_fps = input_pipe(opts.input_fps, optTitle="[-i --input_fps]", verboseIn=opts.verbose)
    
    ######### CUSTOM ############
    if opts.verbose: print "########################################### CUSTOM #########################################"
    #if input_input_fps: print_enumlist(input_input_fps)
    
    
    ### GET INPUT TREE OPTIONS ###
    # INPUT NUM TREES #
    if opts.num_trees: numTbool, numT =  input_int_data(opts.num_trees)
    print numT
    
    # INPUT LEAF NAMES #
    if opts.leaf_names: inLeaves = input_string_data(opts.leaf_names)
    print inLeaves
    
    ### GET NUMBER OF LEAVES ###
    numLeaves = len(inLeaves)
    
    ### CONSTRUCT TREES ###
    tree_random_generator(numTrees=numT, numLeafs=numLeaves, leafIds=inLeaves, outDir=out_path)

    #############################
    if opts.verbose:
        print "############################################ END ###########################################"
        print "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////////////"

##################################################################################################
##################################################################################################
############# GENERAL FUNCTIONS (Required) #######################################################
##################################################################################################
##################################################################################################
#{{{{{{{ TODO IDEAS }}}}}}}#
"""



def input_bool():
    # user input -> check if bool (t\f, T\F, true\false, True\False, TRUE\FALSE, 1\0

Exterior:
- Setup bash_profile so that script runs on 'a'
- Have external bash_profile 'source' script
- Setup private python library with only required dependencies

"""
############################
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

##################################################################################################
##################################################################################################
############# OPTIONAL FUNCTIONS #################################################################
##################################################################################################
##################################################################################################
#################################################
################ STATISTICS #####################


### STATISTICS FUCNTION ### perform OLS regression 
def stat_reg(yIn, xIn):
    regResults = sm.OLS(yIn, xIn).fit()
    print regResults.summary()

#################################################
################## NUMPY ########################!!!

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

#################################################
################# PLOTTING ######################!!!

### PLOT SEABORN ### for dataframe (pandas) - takes in figure axes and plots regression with bootstrapped confidence interval
# orderFit = 1 fits linear regression, greater > 1 fits polynomial regression
# logisticIn = if true assumes y is binary - fits logistic regression 
# robustReg = if true performs robust regression (statsmodels) to de-weight outliers
def plot_reg(xLabel, yLabel, dataFrameIn, 
             axIn=None, nestVar=None, numBootstraps=1000, confidenceInterval=95, orderFit=1, 
             logisticIn=False, robustReg=False, logXIn=False, confoundX=None, confoundY=None):
    stat_reg(dataFrameIn[yLabel],dataFrameIn[xLabel])
    plt.show(sns.regplot(xLabel, yLabel, data=dataFrameIn, ax=axIn, x_ci='ci', scatter=True, fit_reg=True, ci=confidenceInterval, n_boot=numBootstraps, units=nestVar, order=orderFit, logistic=logisticIn, lowess=False, robust=robustReg, logx=logXIn, x_partial=confoundX, y_partial=confoundY, truncate=False, dropna=True, x_jitter=None, y_jitter=None, label=None, color=None, marker='o', scatter_kws=None, line_kws=None))
    

def plot_template():
    ### USER CUSTOMIZATION ###
    # DIMENSIONS #
    xDim = 12
    yDim = 12
    # SUBPLOTS #
    subplotRows=2
    subplotColumns=2
    subplotRowSpacing = 0.2
    subplotColumnSpacing = 0.2
    shareX = False
    shareY = True
    ##########################
    fig, axList = plt.subplots(subplotRows, subplotColumns, figsize=(xDim,yDim), sharex=shareX, sharey=shareY)
    for rowId, rowAx in enumerate(axList):           # LOOP Through Subplot Rows
        for colId, colAx in enumerate(rowAx):        # LOOP Through Subplot Columns
            # INDIVIDUAL PLOT #
            if rowId == 0 and colId == 0:            # EDIT Specific Subplot
                
                ### AXIS CUSTOMIZATION ###
                axisBGColor = 'white'
                xMin = 0
                xMax = 1
                yMin = 0
                yMax = 1
                ###
                colAx.set_axis_bgcolor(axisBGColor)        # Axis background color
                colAx.set_xlim([xMin, xMax])
                colAx.set_ylim([yMin, yMax])
                colAx.set_axis_off()
                
                ##########################
                ### TEXT CUSTOMIZATION ###
                axTextIn = "Title"
                axTextLocXIn = 0.4
                axTextLocYIn = 1.01
                axTextColorIn = 'r'
                axTextSizeIn = 20
                axTextClipIn = False
                axTextRotationIn = 0
                ###
            	colAx = plot_axText(colAx,axText=axTextIn,axTextLocX=axTextLocXIn,axTextLocY=axTextLocYIn,axTextColor=axTextColorIn,axTextSize=axTextSizeIn,axTextClip=axTextClipIn,axTextRotation=axTextRotationIn)
                #colAx = plot_axText(colAx, axText = "Title", axTextLocX = .5, axTextLocY = 1.1, axTextColor = 'k', axTextSize = 20, axTextClip = False, axTextRotation=0)
                ##########################
    
    # DISPLAY & RETURN #
    #fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, hspace=subplotRowSpacing, wspace=subplotColumnSpacing)
    plt.show()
    return

### PLOT FUNCTION ### place text into figure axis and return the axis
def plot_axText(axIn, axText = "Title", axTextLocX = .5, axTextLocY = 1.1, axTextColor = 'k', axTextSize = 20, axTextClip = False, axTextRotation=0):
    axIn.text(axTextLocX,axTextLocY,axText,color=axTextColor,fontsize=axTextSize,clip_on=axTextClip,rotation=axTextRotation)
    return axIn

### MATPLOTLIB ###
def plot_dist( binsIn=10, histIn=True, kdeIn=True, rugIn=True,):
    x = np.random.normal(size=1000)
    sns.distplot(x, bins=binsIn, hist=histIn, kde=kdeIn, rug=rugIn, fit=sp.gamma)
    plt.show()
    return
    
### PLOT SEABORN ### for dataframe (pandas) 
def plot_pairplot(dfIn):
    sns.pairplot(dfIn)
    plt.show()

def plot_jointplot(xLabel, yLabel, dfIn, kindIn="kde", colorIn="m"):
    g = sns.jointplot(x=xLabel, y=yLabel, data=dfIn, kind=kindIn, color=colorIn)
    g.plot_joint(plt.scatter, c="r", s=50, linewidth=1, marker="+")
    g.ax_joint.collections[0].set_alpha(0)
    g.set_axis_labels("$X$", "$Y$")
    plt.show()
    return


def plot_kde():
    x = np.random.normal(size=1000)
    sns.kdeplot(x, shade=True, label="bandwidth: 1")
    sns.kdeplot(x, bw=.2, label="bandwidth: 0.2")
    sns.kdeplot(x, bw=2, label="bandwidth: 2")
    plt.legend();
    
    """
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    plt.xlabel("Component Number", size=20)
    ax.xaxis.set_label_coords(.5, -0.15)

    plt.ylabel("Power", size=20)
    ax.yaxis.set_label_coords(-.09, .5)
    """
#################################################    
################# SEQUENCE ######################!!!
import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.Blast import NCBIWWW


def seq_plot(seqIn, xDim=15,yDim=1):
    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(xDim,yDim), sharex=True)
    print "Sequence Length: "+str(len(seqIn))
    xLen = np.arange(0, len(seqIn))
    pctGC = 0;
    nA = np.zeros(len(seqIn)); nT = np.zeros(len(seqIn)); nC = np.zeros(len(seqIn)); nG = np.zeros(len(seqIn)); nN = np.zeros(len(seqIn))
    for i in xLen:
        ax2.text(i,.75,str(i),color='w',fontsize=8, clip_on=False, rotation=90)
        if seqIn[i] == 'a' or seqIn[i] == 'A': nA[i] = 1; ax1.text(i,.1,"A",color='r',fontsize=16, clip_on=False)
        if seqIn[i] == 't' or seqIn[i] == 'T': nT[i] = 1; ax1.text(i,.1,"T",color='b',fontsize=16, clip_on=False)
        if seqIn[i] == 'c' or seqIn[i] == 'C': nC[i] = 1; ax1.text(i,.1,"C",color='y',fontsize=16, clip_on=False)
        if seqIn[i] == 'g' or seqIn[i] == 'G': nG[i] = 1; ax1.text(i,.1,"G",color='g',fontsize=16, clip_on=False)
        if seqIn[i] == 'n' or seqIn[i] == 'N': nN[i] = 1; ax1.text(i,.1,"N",color='k',fontsize=16, clip_on=False)
	ax1.set_xlim([min(xLen), max(xLen)]); ax1.set_ylim([0, 1]); ax1.set_axis_off()
    ax2.bar(xLen, nA, color='r', label="A")
    ax2.bar(xLen, nT, color='b', label="T")
    ax2.bar(xLen, nC, color='y', label="C")
    ax2.bar(xLen, nG, color='g', label="G")
    ax2.bar(xLen, nN, color='k', label="N")
    ax2.legend(bbox_to_anchor=(1, 1.25), loc=2)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0, wspace=0)
    plt.show()
    
    

# ID: seq
# FUNCTIONS: work with genetic sequence data and files

####################################
 ######## SEQUENCE I/O ########
####################################

### STR_DNA - generate a DNA Sequence ### 
	# Input: String 
	# Return: Sequence 
def str_dna(in_seq):
	global verbose_bool
	dna_seq = Seq(in_seq, IUPAC.unambiguous_dna)
	if verbose_bool is True: print('	* Make Unambiguous DNA Sequence: ' + dna_seq)
	return dna_seq

### STR_RNA - generate a RNA Sequence ### 
	# Input: String
	# Return: Sequence 
def str_rna(in_seq):
	global verbose_bool
	rna_seq = a(in_seq, IUPAC.unambiguous_rna)
	if verbose_bool is True: print('	* Make Unambiguous DNA Sequence: ' + dna_seq)
	return rna_seq

### STR_PROT - generate an amino acid sequence ### Input: String | Return: Sequence 
def str_prot(in_seq):
	global verbose_bool
	prot_seq = Seq(in_seq, IUPAC.protein)
	if verbose_bool is True: print('	* Make Protein Sequence: ' + prot_seq)
	return prot_seq

### STR_NUC - make generic nucleotide sequence ### Input: String | Return: Sequence
def str_nuc(in_seq):
	global verbose_bool
	nuc_seq = Seq(in_seq, generic_nucleotide)
	if verbose_bool is True: print('	* Make Generic Nucleotide Sequence: ' + nuc_seq)
	return nuc_seq
	
### FASTA_SINGLE ### - import a one sequence fasta file
def fasta_single(fasta_fp):
	global verbose_bool
	try:
		record = SeqIO.read(fasta_fp, format="fasta")
		if verbose_bool is True: 
			print('	* Sequence Import from ' + fasta_fp)
		return record
	except NameError:
		print 'Could not import - '+fasta_fp

### FASTA ### - read fasta file into list of sequences
def fasta(fasta_fp):
	global verbose_bool
	list_seqs = []
	if os.path.isfile(fasta_fp):
		for seq_record in SeqIO.parse(fasta_fp, "fasta"):
			if verbose_bool is True: print(seq_record.id), ' - ',(len(seq_record)), ' - ',(repr(seq_record.seq))
			list_seqs.append(seq_record)
	else:
		try:
			for seq_record in SeqIO.parse(os.path.abspath(fasta_fp), "fasta"):
				if verbose_bool is True: print(seq_record.id), ' - ', (len(seq_record)), ' - ', (repr(seq_record.seq))
				list_seqs.append(seq_record)
		except NameError:
			print '	Could Not Import Fasta File! ',fasta_fp
	return list_seqs

### GENBANK ### - read genbank file into list of sequences
# IN: genbank_path
# OUT: dict_sequences
def read_genbank(genbank_fp):
	global verbose_bool
	list_seqs = []
	if os.path.isfile(genbank_fp):
		for seq_record in SeqIO.parse(genbank_fp, "genbank"):
			if verbose_bool is True: print(seq_record.id), ' - ',(len(seq_record)), ' - ',(repr(seq_record.seq))
			list_seqs.append(seq_record)
	else:
		try:
			for seq_record in SeqIO.parse(os.path.abspath(genbank_fp), "genbank"):
				if verbose_bool is True: print(seq_record.id), ' - ',(len(seq_record)), ' - ',(repr(seq_record.seq))
				list_seqs.append(seq_record)
		except NameError:
			print '	Could Not Import genbank File! ',genbank_fp
	return list_seqs

### FASTA_DICT ### - read fasta file into dictionary of sequences
def fasta(fasta_fp):
	global verbose_bool
	dict_seqs = {}
	if os.path.isfile(fasta_fp):
		for seq_record in SeqIO.parse(fasta_fp, "fasta"):
			if verbose_bool is True: print(seq_record.id), ' - ',(len(seq_record)), ' - ',(repr(seq_record.seq))
			dict_seqs[seq_record.id] = seq_record.seq
	else:
		try:
			for seq_record in SeqIO.parse(os.path.abspath(fasta_fp), "fasta"):
				if verbose_bool is True: print(seq_record.id), ' - ', (len(seq_record)), ' - ', (repr(seq_record.seq))
				dict_seqs[seq_record.id] = seq_record.seq
		except NameError:
			print '	Could Not Import Fasta File! ',fasta_fp
	return dict_seqs

### GENBANK_DICT ### - read genbank file into dictionary of sequences
# IN: genbank_path
# OUT: dict_sequences
def read_genbank(genbank_fp):
	global verbose_bool
	dict_seqs = {}
	if os.path.isfile(genbank_fp):
		for seq_record in SeqIO.parse(genbank_fp, "genbank"):
			if verbose_bool is True: print(seq_record.id), ' - ',(len(seq_record)), ' - ',(repr(seq_record.seq))
			dict_seqs[seq_record.id] = seq_record.seq
	else:
		try:
			for seq_record in SeqIO.parse(os.path.abspath(genbank_fp), "genbank"):
				if verbose_bool is True: print(seq_record.id), ' - ',(len(seq_record)), ' - ',(repr(seq_record.seq))
				dict_seqs[seq_record.id] = seq_record.seq
		except NameError:
			print '	Could Not Import genbank File! ',genbank_fp
	return dict_seqs

### FASTA_EDIT ### edit the headers or sequences of a fasta file by looping through
def fasta_edit(fasta_fp, fasta_out):
	f = open(fasta_out, "w")
	if os.path.isfile(fasta_fp):
		# For each sequence one at a time - change headers to incrementing integers and write
		count = 1
		for seq_record in SeqIO.parse(fasta_fp, "fasta"):
			f.write('>'+str(count)+'\n')
			seq = str(seq_record.seq)
			print str(count) + ' -> ' + str(len(seq))
			f.write(seq+'\n')
			count = count + 1
		f.close()


########################################
 ####### SEQUENCE MANIPULATION #######
########################################

### LOOP ### - loop through the characters of a sequence ###
def loop(seq):
	global verbose_bool
	if verbose_bool is True: print('  * Sequence:')
	for index, letter in enumerate(seq):
		if verbose_bool is True: print("        %i %s" % (index, letter))

### LENGTH ### - take in a sequence and return its length ### Input: Sequence | Return: Int
def length(seq):
	global verbose_bool
	if verbose_bool is True: print('  * Sequence Length: '+str(len(seq)))
	return len(seq)

### COUNT_SUBSEQ ### - count how many times a non-overlapping subsequence appears ### Input: Seq, SubSeq
def count_subseq(seq, subseq):
	global verbose_bool
	c = seq.count(subseq)
	if verbose_bool is True: print('  * Subseq ' + subseq + ' appears ' + str(c) + ' times.')
	return c

### COUNT_GC ### - report the GC content of the sequence ### 
def count_gc(seq):
	global verbose_bool
	c = GC(seq)
	if verbose_bool is True: print('  * GC/CG Content ' + str(c) + ' %')
	return c

### CODONS ### - get sequences of each of the three codon positions ###
def codons(seq):
	global verbose_bool
	c1 = seq[0::3]
	c2 = seq[1::3]
	c3 = seq[2::3]
	if verbose_bool is True: print('  * Codon 1 Positions: '+c1)
	if verbose_bool is True: print('  * Codon 2 Positions: '+c2)
	if verbose_bool is True: print('  * Codon 3 Positions: '+c3)
	return c1, c2, c3

### REV ### - return the inverse ###
def rev(seq):
	global verbose_bool
	c = seq[::-1]
	if verbose_bool is True: print('  * Reverse Sequence: ' + c)
	return c

### CONCAT ### - concatenate two sequences ###
def concat(seq1, seq2):
	global verbose_bool
	c = seq1 + seq2
	if verbose_bool is True: print '  * Concatenated Sequence: ' + c
	return c

### UPPER_CASE - return uppercase version of seq ###
def upper_case(seq):
	global verbose_bool
	c = seq.upper()
	if verbose_bool is True: print '  * Upper Case Sequence: ' + c
	return c

### LOWER_CASE - return lowercase version of seq ###
def lower_case(seq):
	global verbose_bool
	c = seq.lower()
	if verbose_bool is True: print '  * Lower Case Sequence: ' + c
	return c

### COMPLEMENT - return the complement of the sequence ###
def complement(seq):
	global verbose_bool
	c = seq.complement()
	if verbose_bool is True: print '  * Complement Sequence: ' + c
	return c
### REV_COMPLEMENT - return the complement of the sequence ###
def rev_complement(seq):
	global verbose_bool
	c = seq.reverse_complement()
	if verbose_bool is True: print '  * Reverse Complement Sequence: ' + c
	return c

### TRANSCRIBE - convert DNA into RNA of Sequence
def transcribe(seq):
	global verbose_bool
	c = seq.transcribe()
	if verbose_bool is True: print '  * Transcribed Sequence: ' + c
	return c

### TRANSCRIBE_RC - convert DNA into RNA of Reverse Complement Sequence 
def transcribe_rc(seq):
	global verbose_bool
	c = seq.reverse_complement().transcribe()
	if verbose_bool is True: print '  * Transcribed Sequence of Reverse Complement: ' + c
	return c

### TRANSCRIBE_REV - convert RNA back into DNA of Sequence 
def transcribe_rev(seq):
	global verbose_bool
	c = seq.back_transcribe()
	if verbose_bool is True: print '  * Reverse Transcribed Sequence: ' + c
	return c

### TRANSLATE - convert DNA or RNA to Amino Acid Sequence
def translate(seq):
	global verbose_bool
	c = seq.translate()
	if verbose_bool is True: print '  * Translated Sequence: ' + c
	return c

### BLAST ### fasta format -> result handle
def blast(seq):
	try:
		result_handle = NCBIWWW.qblast("blastn", "nt", seq.format("fasta"))
		return result_handle
	except NameError:
		print 'Could not blast and save '+i

### BLAST_SAVE_RAW ### result_handle -> xml file
def blast_save_raw(result_handle, save_loc):
	try:
		save_file = open(save_loc, "w")
		in_d = result_handle.read()
		save_file.write(in_d)
		save_file.close()
		result_handle.close()
		print 'Results of BLAST saved to - '+i+"_blast.xml"	
	except NameError:
		print 'Could not save BLAST results to '+save_loc

#################################################
################## TREES ########################
### RANDOM TREE TOPOLOGY GENERATOR ### 
from ete2 import Tree
def tree_random_generator(numTrees=100, numLeafs=6, leafIds=['A','B','C','D','E','F'], outDir=''):
    for i in np.arange(numTrees):
        t = Tree()
        leafIds = randomize(leafIds)
        t.populate(numLeafs, names_library=leafIds)
        print outDir+'tree_'+str(numLeafs)+'_'+str(i)+'.newick'
        #print t
        t.write(outfile=outDir+'tree_'+str(numLeafs)+'_'+str(i)+'.newick')
    
# Randomize Input List Function
import random
def randomize(a):
    b = []
    for i in range(len(a)):
        element = random.choice(a)
        a.remove(element)
        b.append(element)
    return b

#################################################
################## CUSTOM #######################
#################################################


#################################################
#################### END ########################
#################################################
       
if __name__ == "__main__":
    main()