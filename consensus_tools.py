## CONSENSUS TOOLS
# Description: A set of tools for reconciling multiple variants of a string of text using sequence alignment and consensus algorithms
# Author: Junying Lim
# Date: 13th Jan 2014

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
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Align import AlignInfo

mafft = "/usr/local/bin/mafft"

def create_variant_dict(accession, field, data):
    ''' Creates a dictionary from different transcriptions

        Arguments:
        accession -- string, define unique unique ID field
        field     -- string, define target field to build a dictionary of
        data      -- pandas.core.frame.Dataframe object, 
                     must contain specified accession and field columns

        Returns:
            defaultdict object
    '''

    # Data checks
    if not isinstance(data, pd.DataFrame):
        raise Exception("Data must be a pandas.core.frame.DataFrame object")
    if accession not in list(data.columns):
        raise Exception("Accession field not found in data object")
    if field not in list(data.columns):
        raise Exception("Target field not found in data object")

    # Generate dictionary of accessions paired with list of field entries
    entry_list = zip(list(data[accession]), list(data[field]))
    entry_id = defaultdict(list)
    for k,v in entry_list:
        entry_id[k].append(v)

    return entry_id

def best_transcript(accession,  field, data, method):
    ''' Selects the best variant based on its similarity to other variants
        
        Arguments:
        accession   -- string, define unique ID field
        field       -- string, define target field to resolve,
        data        -- pandas.core.frame.Dataframe object,
                    must contain specified accession and field columns
        method      -- either fuzzy string matching ('fuzzy')
                    or Levenshtein ('distance') method,

        Returns:
        a pandas.core.frame.Dataframe object

    '''

    # Data checks
    if method not in ["fuzzy", "distance"]:
        raise Exception("Method argument not recognized. Must be either 'fuzzy' or 'distance'")

    # Generate dictionary of accessions paired with list of field entries
    entry_id = create_variant_dict(accession, field, data)

    # Reconcile entries
    if method == "fuzzy":
        # Create a new dictionary to hold results
        entry_fuzzy_results = defaultdict(list)
        entry_pairs = defaultdict(list)
        for k, v in entry_id.items():
            pairs = []
            for pair in it.combinations(v, 2):
                # pair = map(str.lower, pair)
                entry_pairs[k].append(pair)
                
                # Calculate the "fuzz" score
                c=int(fuzz.token_sort_ratio(str(pair[0]),str(pair[1])))
                entry_fuzzy_results[k].append(c)

        acc_id = []
        entry_est = []
        for acc, c_scores in entry_fuzzy_results.items():
            max_c = c_scores.index(max(c_scores))
            estimate = entry_pairs[acc][max_c][1]
            estimate = estimate.capitalize()
            acc_id.append(acc)
            entry_est.append(estimate)
                        
        return pd.DataFrame({str(accession): acc_id, str(field): entry_est})

    elif method == "distance":
        # Create a new dictionary to hold results
        entry_dist_results = defaultdict(list)
        entry_pairs = defaultdict(list)
        # For each specimen:
        for k, v in entry_id.items():
            pairs = []
            for pair in it.combinations(v, 2):
                # pair = map(str.lower, pair)
                entry_pairs[k].append(pair)
                
                if str(pair[0]) == "" and str(pair[1]) == "":
                    L_dist = int(100) # Arbitrarily large number

                else:
                    L_dist = int(distance(str(pair[0]), str(pair[1])))

                entry_dist_results[k].append(L_dist)
                         
        acc_id = []
        entry_est = []
        for acc, L_scores in entry_dist_results.items():
            min_L = L_scores.index(min(L_scores))
            estimate = entry_pairs[acc][min_L][1]
            acc_id.append(acc)
            entry_est.append(estimate)
                        
        return pd.DataFrame({str(accession): acc_id, str(field): entry_est})

def token_align(x, wdir, consensus_method):
    ''' Implements a sequence alignment and consensus finding on tokenized strings

        Arguments:
        x       -- List of strings

        Returns:
        Consensus string

    '''
    if consensus_method not in ["dumb", "dumber"]:
        raise Exception("Consensus method not recognized. Must be either 'dumb' or 'dumber'")

    
    # Concatenate strings and tokenize 
    regexptoken = nltk.RegexpTokenizer(pattern = '\s+', gaps = True)                           
    y = string.join(x)
    tokens = regexptoken.tokenize(y)
    # tokens = nltk.word_tokenize(y) #word_tokenize assumes you are a sentence, so any ultimate periods are tokenize separately

    # Find unique tokens
    unique_tokens = []
    for token in tokens:
        if token not in unique_tokens:
            unique_tokens.append(token)
            
    # Create a dictionary of tokens
    unique_token_id = dict()
    index = 0
    for unique_token in unique_tokens:
        unique_token_id[unique_token] = string.ascii_letters[index]
        index += 1

    # Convert strings into token IDs
    entry_token_strings = []
    for entry in x:
        entry_tokens = regexptoken.tokenize(entry) #nltk.word_tokenize(entry)
        entry_token_string = ""
        for token in entry_tokens:
            #print token
            entry_token_string += unique_token_id[token]
            #print entry_token_string
        entry_token_strings.append(entry_token_string)

    # Convert strings into SeqRecord objects
    entry_token_strings = [Seq(token) for token in entry_token_strings]
    temp = [SeqRecord(entry_token_string, id = str(acc)) for acc, entry_token_string in enumerate(entry_token_strings)]

    # Write a fasta file from the list of 
    temp_file = os.path.join(wdir, "temp.fasta")
    SeqIO.write(temp, temp_file, "fasta")

    # Align using alignment algorithm MAFFT
    res = subprocess.check_output([mafft, '--text','--localpair','--maxiterate','1000',temp_file])
    res = res.decode("utf-8")

    # Export results into working dir
    out_file = os.path.join(wdir, "temp_align.fasta")  # create alignment file name

    f = open(out_file, 'wb')
    f.write(bytes(res, "UTF-8"))
    f.close()

    # Import alignment
    alignres = AlignIO.read(out_file, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignres)

    # Determine consensus
    consensus = alignment_consensus(summary_align, method = consensus_method)

    # Reinterpret consensus
    inv_unique_token_id = {v:k for k, v in unique_token_id.items()}
    consensus = [inv_unique_token_id[token] for token in consensus]
    consensus = string.join(consensus, sep = " ")

    # Remove temporary file
    os.remove(out_file)
    os.remove(temp_file)
    
    return str(consensus)


def character_align(x, wdir, consensus_method):
    ''' x   -- list of strings
    
        Returns: Single string
    '''
    if consensus_method not in ["dumb", "dumber"]:
        raise Exception("Consensus method not recognized. Must be either 'dumb' or 'dumber'")

    
    # Placeholder characters (= * < > ( ) not allowed in mafft)
    y = [re.sub("\s", "_", string) for string in x] # convert spaces into underscores
    y = [re.sub("=", "", string) for string in y]
    y = [re.sub("<", "", string) for string in y]
    y = [re.sub(">", "", string) for string in y]
    y = [re.sub("\*", "", string) for string in y]
    y = [re.sub("\.", "%", string) for string in y]
    y = [re.sub("\(", "[", string) for string in y]
    y = [re.sub("\)", "]", string) for string in y]
    
    # Convert string into SeqRecord objects
    temp = [Seq(string, "alphabet") for string in y]
    temp = [SeqRecord(string, id = str(acc)) for acc, string in enumerate(temp)]
    temp_file = os.path.join(wdir, "temp.fasta")
    SeqIO.write(temp, temp_file, "fasta")

    # Subprocessing MAFFT
    res = subprocess.check_output([mafft, '--text','--localpair','--maxiterate','1000', temp_file])
    res = res.decode("utf-8")

    # Export results into working dir
    out_file = os.path.join(wdir, "temp_align.fasta")  # create alignment file name

    f = open(out_file, 'wb')
    f.write(bytes(res, "UTF-8"))
    f.close()

    # Import alignment
    alignres = AlignIO.read(out_file, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignres)

    # Determine consensus
    consensus = alignment_consensus(summary_align, method = consensus_method)

    # Strips excessive white space
    consensus = consensus.strip()

    # Back convert placeholders
    consensus = re.sub("_", " ", str(consensus))
    consensus = re.sub("%", ".", str(consensus))
    consensus = re.sub("\[", "(", str(consensus))
    consensus = re.sub("\]", ")", str(consensus))

    # Remove temporary file
    os.remove(out_file)
    os.remove(temp_file)

    # Return
    return consensus

def alignment_consensus(alignment, method):

    '''
    Takes a summary alignment object

    Args:

    '''
    if method not in ["dumb", "dumber"]:
        raise Exception("Consensus method not recognized. Must be either 'dumb' or 'dumber'")
    
    if method == "dumb":
        consensus = alignment.dumb_consensus(threshold = 0.5, require_multiple = 1, consensus_alpha = None, ambiguous = "")
    elif method == "dumber":
        consensus = dumber_consensus(alignment, threshold = 0.5, ambiguous = "")

    return consensus

def variant_consensus(accession, field, data, align_method, consensus_method, wdir):
    if align_method not in ["character", "token"]:
        raise Exception("Alignment method not recognized. Must be either 'fuzzy' or 'distance'")

    if consensus_method not in ["dumb", "dumber"]:
        raise Exception("Consensus method not recognized. Must be either 'dumb' or 'dumber'")

    
    entry_id = create_variant_dict(accession, field, data)

    print("\nImplementing consensus procedure on", field, "field, using", align_method, "alignment method and", consensus_method, "consensus method")

    # Find consensus in NfN data
    entry_results = defaultdict(list)
    for k,v in entry_id.items():
        print("Reconciling transcriptions for", k)

        # If entries are identical, then entry is consensus
        if len(set(v)) == 1:
            entry_results[k].append(v[0])

        # If all entries are one character in length or below, then consensus is probably nothing
        elif sum([len(i) < 2 for i in v]) == len(v):
            entry_results[k].append("")

        # If entries are not identical, use consensus
        else:
            if align_method == "character":
                entry_results[k].append(character_align(v, wdir, consensus_method))
            elif align_method == "token":
                entry_results[k].append(token_align(v, wdir, consensus_method))

    # Convert results into dataframe
    est = [str(est[0]) for est in entry_results.values()] # Necessary to index 0 and default dict values are lists
    acc = [str(acc) for acc in entry_results.keys()]
    results = pd.DataFrame({str(accession):acc, str(field):est})

    # Export
    return results

def vote_count(accession, field, data):

    ''' Uses a vote counting procedure for selecting the best. For fields with discrete states.
        
        Arguments:
        accession   -- string, define unique ID field
        field       -- string, define target field to resolve,
        data        -- pandas.core.frame.Dataframe object,
                    must contain specified accession and field columns

        Returns:
        a pandas.core.frame.Dataframe object
    '''
    print("Implementing vote-counting procedure on", field, "field.")
    entry_id = create_variant_dict(accession, field, data)

    entry_results = defaultdict(list)
    for k,v in entry_id.items():
        print("Reconciling transcriptions for", k)
        count_vote = [v.count(i) for i in set(v)]
        max_vote = list(set(v))[count_vote.index(max(count_vote))]
        entry_results[k].append(max_vote)

    vote = [str(v[0]) for v in entry_results.values()]
    key = [str(k) for k in entry_results.keys()]
    results = pd.DataFrame({str(accession):key, str(field):vote})
    return results


def dumber_consensus(self, threshold, ambiguous = ""):
        # find the length of the consensus we are creating
        con_len = self.alignment.get_alignment_length()
        num_atoms = len(self.alignment)
        consensus = ""
        # go through each seq item
        for n in range(con_len):
            # keep track of the counts of the different atoms we get
            atom_dict = {}
            
            for record in self.alignment._records:
                # make sure we haven't run past the end of any sequences
                # if they are of different lengths
                if n < len(record.seq):
                    if record.seq[n] != '-' and record.seq[n] != '.':
                        if record.seq[n] not in atom_dict.keys():
                            atom_dict[record.seq[n]] = 1
                        else:
                            atom_dict[record.seq[n]] = \
                              atom_dict[record.seq[n]] + 1

            max_atoms = []
            max_size = 0
            for atom in atom_dict.keys():
                if atom_dict[atom] > max_size:
                    max_atoms = [atom]
                    max_size = atom_dict[atom]
                elif atom_dict[atom] == max_size:
                    max_atoms.append(atom)
                    
            # max_atoms now = most frequently occuring character
            # max_size now = no. of occurences of the charater
            # atom_dict = carries all variants

            # If there is a clear conesnsus character, append that character 
            if (len(max_atoms) == 1) and ((float(max_size)/float(num_atoms))
                                         >= threshold):
                consensus = consensus + max_atoms[0]
    
            # If two characters are tied and both are characters, append the lower-case character (usually equally good, e.g., lowercase and uppercase) 
            elif len(max_atoms) == 2 and sum([max_atom in string.ascii_letters for max_atom in max_atoms]) == 2:
                if max_atoms[0] in string.ascii_lowercase:
                    consensus = consensus + max_atoms[0]
                elif max_atoms[1] in string.ascii_lowercase:
                    consensus = consensus + max_atoms[1]

            # Else, append an ambiguous character
            else:
                consensus = consensus + ambiguous

        return consensus


def metadata_handling(accession, field, data, delim = "|"):

    ''' Summarizes metadata fields by simply appending and delimiting them into a string
        
        Arguments:
        accession   -- string, define unique ID field
        field       -- string, define target field to resolve,
        data        -- pandas.core.frame.Dataframe object,
                    must contain specified accession and field columns

        Returns:
        a pandas.core.frame.Dataframe object
    '''
    
    # Create a dictionary of accession keys to metadata values
    entry_id = create_variant_dict(accession, field, data)
    entry_results = defaultdict(list)

    for k,v in entry_id.items():
        entry_results[k].append(delim.join(v)) #[j for j in set(v)[0]] this was for picking just 1
   
    key = [str(k) for k in entry_results.keys()]
    meta = [str(v[0]) for v in entry_results.values()]
    results = pd.DataFrame({str(accession):key, str(field):meta})
    return results


def reflist_check(datalst, reflst, threshold):

    ''' Uses fuzzy string matching to identify the most likely name from a reference list
        
        Arguments:
        datalst     -- list, contains names to check
        reflst      -- list, contains list of possible names

        Returns:
        a tuple containing a list of best matches, and certainty values

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


def refcheck(datalst, reflst, threshold):
    estimate = list()
    score = list()

    
    counter = 0
    for x in datalst:           
        if counter % 100 == 0:
                print(counter, "of", len(datalst), "entries normalized.")

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





