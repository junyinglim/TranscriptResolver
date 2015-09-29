
import csv
import os
import numpy as np
import re # regular expressions
import argparse #For command line arguments
import pandas as pd # data frame functionality
import pymysql

## MAIN ##
def main():
    args = parser.parse_args()
    
    ## Preamble ========================
    print("\n\n\n")
    print("=" * 50)
    print("PREPARING TRANSCRIPTS FOR RESOLUTION!!")
    
    ## Import transcription file ========================
    print("\nImporting transcription file ", args.file, " ...")
    data = pd.read_csv(args.file, encoding = "ISO-8859-1", dtype = "object")
    data = data.fillna("")

    ## Prepping transcription file ========================
    print("\nPrepping transcription file...")

    # Exclude non-calbug entries
    print("\nExcluding non-Calbug transcriptions ...")
    data = data[data["collection"] == "Calbug"] 

    # Convert Collector names to lower case to facilitate alignments
    data["Collector"] = [collector.lower() for collector in data["Collector"]]
    
    # Connect to essig database
    conn = pymysql.connect(host = "gall.bnhm.berkeley.edu",\
                           user = args.username, # should be args.username
                           passwd = args.password, # should be args.password
                           db = "essig")
    essigIDs = pd.read_sql('select bnhm_id from eme;', con=conn)
    essigIDs = list(essigIDs["bnhm_id"])
    
    
    # exclude specimens that are already in the database
    ##comment## are all the nfn entries meant for calbug, this stage might basically exclude all the non-essig transcirptions
    ##comment## assumes that a specimen ID does not have multiple filenames
    
    print("\nExcluding specimens that have already been databased ...")
    unique_filename = list(set(data[args.col_id])) # create a dictionary of bnhm_id and filename
    id_dict = dict()

    for f in unique_filename:
        ids = re.sub("\s+", "", f)
        ids = re.search("((EMEC|LACMENT|CASENT|UCBME|CIS|UCRCENT|SDNHM|UMMZI|SBMNHENT)[0-9]*)", ids).group()
        id_dict[ids] = f
    
    unique_ids = set(id_dict.keys())
    completed_ids = list(unique_ids.intersection(essigIDs))

    completed_filenames = [id_dict[f] for f in completed_ids]
    completed_filenames

    # exclude filenames whose bnhm_id is already in the essig database
    data = data[~data["filename"].isin(completed_filenames)] # ~ means the inverse; so filenames that have not been completed
    
    ## Export the new data file ========================
    if args.output:
        outputfile = args.output
    else:
        outputfile = "transcript_prepped"
    
    print("\nExporting prepared transcriptions to", os.getcwd())
    data.to_csv(outputfile + ".csv", index = False)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="transcriptResolver - Let's resolve some crowd-sourced transcripts!")
    parser.add_argument("-file", "-f", help = "File with transcriptions. Full file path required!")
    parser.add_argument("-output", help = "Output file name. Full path not required")
    parser.add_argument("-col_id", help = "List of columns to be checked")
    parser.add_argument("-username", help = "Username. Access to essig SQL database")
    parser.add_argument("-password", help = "Password. Access to essig SQL database")
    main()
    
##todo## logging the results