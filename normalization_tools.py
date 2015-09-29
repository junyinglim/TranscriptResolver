## NORMALIZATION TOOLS
# Description: A set of tools for fuzzy matching crowd-sourced verbatim to refernece libraries
# Author: Junying Lim
# Date: 21th Sept 2015

# Notes:
# Developed during the CITScribe hackathon organized by iDigBio, Gainsville, Florida (15 Dec - 20 Dec 2013)
# For more information on using MAFFT for string matching http://mafft.cbrc.jp/alignment/software/textcomparison.html


## DEPENDENCIES
import pandas as pd 
from collections import defaultdict
import itertools as it
from fuzzywuzzy import fuzz # Fuzzy string matching for best_transcript()
from fuzzywuzzy import process # For reflst matching
from Levenshtein import * # Levenshtein distance for best_transcript() 

import os # Path tools
import string # String tools
import nltk # For tokenizing
import subprocess # For subprocessing MAFFT
import re # Regular expressions


def fill_pd(x, pd):
    ''' Populates a pandas.DataFrame by row

        Args:
            x : list of lists (each element of those lists)
            pd: empty pandas data frame

        Returns:
            A pandas.DataFrame

    ''' 
    row = 0
    for i in x:
        col = 0

        if len(i) > len(pd.columns):
            print("There were more entries than columns, will leave row blank")
            pd.loc[row, col] = ""

        # Else fill it up as per normal
        else:
            for j in i:
                pd.loc[row, list(pd.columns.values)[col]] = str(j)
                col += 1

        row += 1

    return pd.fillna("")


def refcheck(datalst, reflst, threshold):
    ''' Uses fuzzy string matching to identify the most likely name from a reference list
        
        Args:
            datalst     : list, contains names to check
            reflst      : list, contains list of possible names
            threshold   : threshold similarity before a name on reference list is considered likely. If 1, then matches must be exact.

        Returns:
            a tuple of 2 lists - one list of the best matches, and the second of certainty values

    '''
    estimate = list()
    score = list()

    counter = 0
    for x in datalst:           
        #if counter % 100 == 0:
         #       print(counter, "of", len(datalst), "entries normalized.")

    # If empty, then return empty
        if x == "":
            estimate.append("")
            score.append("NA")
            counter += 1

    # Else, find the most similar collector name from the whole reference list 
        else:
            temp = [ratio(x, y) for y in reflst]
            max_ratio = max(temp)
                        
            if max_ratio < threshold:        
                estimate.append("")
                score.append(max_ratio)

            else:               
                estimate.append(reflst[temp.index(max_ratio)])
                score.append(max_ratio)

            counter += 1

    return (estimate, score)







def reflist_check(datalst, reflst, threshold):

    ''' Uses fuzzy string matching to identify the most likely name from a reference list
        
        Args:
            datalst     : list, contains names to check
            reflst      : list, contains list of possible names

        Returns:
            a tuple containing a list of best matches, and certainty values

    >>> x = ["E.G.Lyndsey"]
    >>> y = ["Gordon", "Lyndsey", "E.G.Lyndsey"]
    >>> reflist_check(x,y)

    '''
    # Create empty lists
    estimate = ["NA"] * len(datalst)
    certainty = ["NA"] * len(datalst)
    reflistkeys = [key.decode('latin-1') for key in reflst.keys()]

    counter = 0
    for x in datalst:

        # Print counter (the function takes a long time)
        if counter % 100 == 0:
            print(counter, "of", len(datalst), "entries normalized.")

        if x == "":
            estimate[counter] = ""
            certainty[counter] = "NA"
            counter += 1

        else:
            
            temp = process.extractOne(x, reflistkeys)
            
            if temp[1] < threshold:
                estimate[counter] = ""
                certainty[counter] = temp[1]

            else:
                estimate[counter] = reflst[temp[0]]
                certainty[counter] = temp[1]
            counter += 1

    return (estimate, certainty)
