
import csv
import os
import itertools as it
import numpy as np
import re # regular expressions
import argparse #For command line arguments

from consensus_tools import * # custom functions to run transcript resolving
from collections import defaultdict # utility functions to create dictionaries
import pandas as pd # data frame functionality
from fuzzywuzzy import process, fuzz # Functions that are useful for fuzzy string matching (https://github.com/seatgeek/fuzzywuzzy)
from functools import reduce # for the reduce function

import webbrowser

class transcriptResolver:
    def __init__(self, args): # __init__ always run when an instance of the class is created

        ## Define stem name ========================
        if args.stem:
            self.stem = args.stem
        else:
            print("\nPlease input a 'stem' name to act as a prefix to all output (e.g., stemname_calbug.csv)")
            self.stem = input("Stem name: ")
        
        print("\nStem name for all outputs will be '" + self.stem + "'")
        
        ## Define working directory ========================
        if args.wd:
            self.wd = args.wd
        else:
            temp = input("Working directory: ")
            self.wd = temp

        print("\nUsing working directory '" + self.wd + "' ...")

        ## Define transcription file ========================
        if args.file:
            tempfile = args.file
        else:
            tempfile = input("\nInput your working file name. \nWorking file should be in your stated working directory:")
        
        filedir = os.path.join(self.wd, tempfile)
        print("\nFile directory will be '" + filedir + "'")

        
        ## Define id column ========================
        ##todo## need to check if ID column is in the file
        if args.col_id:
            self.col_id = args.col_id
        else:
            self.col_id = input("\nDefine the column name specifying unique IDs (e.g., specimen ID):")

        print("\nColumn name that specifies unique ID is " + self.col_id)
        
        ## Define target columns ========================
        ##todo## need to check if target columns are in the file
        if args.col_target and args.col_method:
            temp_target = args.col_target.strip("[|]").split(",")
            temp_method = args.col_method.strip("[|]").split(",")
            if len(temp_target) > 0 and len(temp_method) > 0 and len(temp_method) == len(temp_target):
                self.col_target = temp_target
                self.col_method = temp_method
            
        else:
            while True:
                temp_target = input("\nDefine the column name to be resolved:\n (Enter nothing to continue to the next step)")
                temp_method = input("\nPlease define the method for which you would like to use on this column: \n(Enter nothing to continue to the next step)")
            
                if(temp_method == "" | temp_target == ""):
                    break
            
                else:
                    self.col_target = temp_target
                    self.col_method = temp_target
        
        [print("\nUsing method", y, "for column", x) for x, y in zip(self.col_target, self.col_method)]
        
        ## Import file ========================
        self.file = pd.read_csv(filedir, dtype=object,\
                                usecols = self.col_target + self.col_id) # only use columns that were supplied
        
        ##todo## might be possible problems with encoding? if yes, then supply encoding = "ISO-8859-1" to read_csv()
        self.file = self.file.fillna("") # Converts all NaNs into empty strings for alignment
        
        
## MAIN ##
def main():
    args = parser.parse_args()
    if args.version:
        print("v1.0")
    elif args.manual:
        webbrowser.open("https://github.com/junyinglim/Notes-from-Nature")
    else:
        print("\n\n\n")
        print("=" * 50)
        print("WELCOME TO TRANSCRIPT RESOLVER!!")
        print("Let's resolve some replicate transcripts! \n")
        print("Please visit https://github.com/junyinglim/Notes-from-Nature \nfor a short explanation of the transcript resolution methods available \n")
        print("This crude program was written by Jun Ying Lim (junyinglim@gmail.com) \nfor the Essig Museum of Entomology at UC Berkeley")
        print("=" * 50)
        print("\n\n\n")
               
    ## Startup
    
    currentArgs = transcriptResolver(args)
    print(currentArgs.file)
    
    # Create empty list
    results = []
    for col_no in range(len(currentArgs.col_target)):
        if currentArgs.col_method[col_no] == "vote_count":
            df = vote_count(accession = currentArgs.col_id,\
                            field = currentArgs.col_target[col_no],\
                            data = currentArgs.file)
        elif currentArgs.col_method[col_no] == "consensus":
            df = variant_consensus(accession = currentArgs.col_id,\
                                   field = currentArgs.col_target[col_no],\
                                   align_method = "character",\
                                   consensus_method = "dumber",\
                                   wdir = currentArgs.wd,\
                                   data = currentArgs.file)
        elif currentArgs.col_method[col_no] == "metadata":
            df = metadata_handling(accession = currentArgs.col_id,\
                                   field = currentArgs.col_target[col_no],\
                                   data = currentArgs.file)
        else:
            ##todo## write a proper error handling here
            print("Sorry, method supplied is not valid")
        
        # Add data frame to the results list
        results.append(df)
        
    # Merge results
    allResults = reduce(lambda a, d: pd.merge(a, d, on = currentArgs.col_id), results)
            
    finalDir = os.path.join(currentArgs.wd, currentArgs.stem + "_results.csv")
    allResults.to_csv(finalDir, index = False)
    print("\nExporting results to", finalDir)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="phyloGenerator - phylogeny generation for ecologists.", epilog="Help at http://willpearse.github.com/phyloGenerator - written by Will Pearse")
    parser.add_argument("--version", action="store_true", help="Display version information.")
    parser.add_argument("--manual", action="store_true", help="(Attempt to) open browser and show help")
    parser.add_argument("-stem", "-n", help="'Stem' name for all output files.") # for command line
    parser.add_argument("-wd", help = "Working directory")
    parser.add_argument("-file", "-f", help = "File with transcriptions")
    parser.add_argument("-col_id", help = "List of columns to be resolved")
    parser.add_argument("-col_target", help = "Target column. Must be in the format -col_target [target1,target2,target3]")
    parser.add_argument("-col_method", help = "Method. Must be in the format -col_method [method1,method2,method3]")
    main()
    
##todo## logging the results