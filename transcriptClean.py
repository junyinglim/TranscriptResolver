
import os
import pandas as pd
import re
import numpy as np
import itertools as it
import csv

from name_splitter import * # Code courtesy of Charles McCallum
from normalization_tools import *
import pymysql
import argparse
from functools import reduce # for the reduce function


class transcriptCleaner:
    def __init__(self, args):
        
        ## IMPORTING FILE ========================
        print("\nImporting file")
        if args.wd:
            data = pd.read_csv(os.path.join(args.wd, args.file), encoding = "ISO-8859-1", dtype = "object")
        else:
            data = pd.read_csv(args.file, encoding = "ISO-8859-1", dtype = "object")
        
        self.data = data.fillna("")
        
        ## CREATING REFERNECE LISTS ========================
        print("\nCreating reference lists for some fields")
        
        conn = pymysql.connect(host = "gall.bnhm.berkeley.edu",
                               user = args.username, # should be args.username
                               passwd = args.password, # should be args.password
                               db = "essig")

        self.essig_collector     = pd.read_sql("select name_full, name_short, collector from eme_people where collector = 1", con = conn)
        
        #essig_canprov       = pd.read_sql("select * from canadian_provinces", con = conn)
        #essig_country       = pd.read_sql("select * from country", con = conn)
        #essig_mexstate      = pd.read_sql("select * from mexican_states", con = conn)
        #essig_state         = pd.read_sql("select * from state", con = conn)
        #essig_county        = pd.read_sql("select * from county", con = conn)
        
        self.essig_canprov = pd.read_csv(os.path.join(os.getcwd(), "reference", "essig_canprov.csv"))
        self.essig_country = pd.read_csv(os.path.join(os.getcwd(), "reference", "essig_country.csv"), encoding = "latin-1")
        self.essig_mexstate = pd.read_csv(os.path.join(os.getcwd(), "reference", "essig_mexstate.csv"), encoding = "latin-1")
        self.essig_statecounty = pd.read_csv(os.path.join(os.getcwd(), "reference", "essig_statecounty.csv"))
        self.essig_holdinginst   = pd.read_csv(os.path.join(os.getcwd(), "reference", "essig_inst.csv"))
    
    def normalizeCollector(self):
        print("\nSplitting collector name strings ...")
    
        # Create empty data frame
        columns = ["Collector_split", "Collector2_split", "Collector3_split", "Collector4_split", "Collector5_split"] 
        index = range(len(self.data))
        split_collector = pd.DataFrame(index=index, columns=columns)

        # Split names
        namesplit = NameSplitter()
        split_collectors =  [namesplit.split(collectors) for collectors in list(self.data['Collector'])]

        # Fill new data frame
        split_collector = fill_pd(split_collectors, split_collector)

        # NORMALIZE COLLECTOR NAMES ========================    
        print("\nNormalizing collector names based on reference list")

        # Build hash table of collectors
        self.essig_collector = self.essig_collector.fillna("")
        self.essig_collector['name_short'] = [col.strip("=") for col in list(self.essig_collector['name_short'])] # Remove the "=" characters in the essig database

        essig_ref = dict()
        for full_name, name_short in zip(self.essig_collector['name_full'], self.essig_collector['name_short']):
            #full_name = full_name.encode("latin-1")
            #name_short = name_short.encode("latin-1")
            
            if full_name != "":
                essig_ref[re.sub('[\W_]+', '', full_name).lower()] = full_name # Remove spaces and non alphanumeric characters from reference list

            if name_short != "":
                temp = name_short.split(", ")
                for sec_name in temp:                    
                    essig_ref[re.sub('[\W_]+', '', sec_name).lower()] = full_name # Remove spaces and non alphanumeric characters from reference list

        # To deal with non-ascii characters
        essig_ref_keys = list(essig_ref.keys())

        # Remove spaces and non alphanumeric characters from collector names
        split_collector['Collector_split'] = [re.sub('[\W_]+', '', col).lower() for col in split_collector['Collector_split']]
        split_collector['Collector2_split'] = [re.sub('[\W_]+', '', col).lower() for col in split_collector['Collector2_split']]
        split_collector['Collector3_split'] = [re.sub('[\W_]+', '', col).lower() for col in split_collector['Collector3_split']]
        split_collector['Collector4_split'] = [re.sub('[\W_]+', '', col).lower() for col in split_collector['Collector4_split']]
        split_collector['Collector5_split'] = [re.sub('[\W_]+', '', col).lower() for col in split_collector['Collector5_split']]

        # Check collector names
        threshold = 0.9 # Specify similarity threshold
        clean_collector1 = refcheck(split_collector['Collector_split'], essig_ref_keys, threshold = threshold)
        clean_collector2 = refcheck(split_collector['Collector2_split'], essig_ref_keys, threshold = threshold)
        clean_collector3 = refcheck(split_collector['Collector3_split'], essig_ref_keys, threshold = threshold)
        clean_collector4 = refcheck(split_collector['Collector4_split'], essig_ref_keys, threshold = threshold)
        clean_collector5 = refcheck(split_collector['Collector5_split'], essig_ref_keys, threshold = threshold)

        # Create new collector data frame
        index = range(len(self.data))
        clean_collector = pd.DataFrame(index = index,\
                                       columns = ["Collector", "Collector2", "Collector3", "Collector4", "Collector5"])
        clean_collector_cert = pd.DataFrame(index = index,\
                                            columns = ["Collector_cert", "Collector2_cert", "Collector3_cert", "Collector4_cert", "Collector5_cert"])

        # Output normalized collectors and certainty levels
        clean_collector["Collector"] = [essig_ref[key] if key is not "" else "" for key in clean_collector1[0]]
        clean_collector["Collector2"] = [essig_ref[key] if key is not "" else "" for key in clean_collector2[0]]
        clean_collector["Collector3"] = [essig_ref[key] if key is not "" else "" for key in clean_collector3[0]]
        clean_collector["Collector4"] = [essig_ref[key] if key is not "" else "" for key in clean_collector4[0]]
        clean_collector["Collector5"] = [essig_ref[key] if key is not "" else "" for key in clean_collector5[0]]

        clean_collector_cert["Collector_cert"] = clean_collector1[1]
        clean_collector_cert["Collector2_cert"] = clean_collector2[1]
        clean_collector_cert["Collector3_cert"] = clean_collector3[1]
        clean_collector_cert["Collector4_cert"] = clean_collector4[1]
        clean_collector_cert["Collector5_cert"] = clean_collector5[1]

        # Export collectors
        self.clean_collector = clean_collector
        self.clean_collector_cert = clean_collector_cert

    def normalizeDates(self):
        index = range(len(self.data))

        # Create empty date dataframes
        split_begin_date = pd.DataFrame(index = index, columns = ["MonthCollected", "DayCollected" , "YearCollected"] )
        split_end_date = pd.DataFrame(index = index, columns = ["MonthCollected2","DayCollected2", "YearCollected2"])

        # String split dates and populate dataframes
        split_begin_dates = [date.split("/") for date in self.data["Begin Date Collected"]]
        split_begin_date = fill_pd(x = split_begin_dates, pd = split_begin_date)
        split_end_dates = [date.split("/") for date in self.data["End Date Collected"]]
        split_end_date = fill_pd(x = split_end_dates, pd = split_end_date)

        # Deal with zeros (should coerced to be an empty field)
        split_begin_date["MonthCollected"] = [date.replace('00', '') for date in split_begin_date["MonthCollected"]]
        split_begin_date["DayCollected"] = [date.replace('00', '') for date in split_begin_date["DayCollected"]]
        split_begin_date["YearCollected"] = [date.replace('0000', '') for date in split_begin_date["YearCollected"]]

        split_end_date["MonthCollected2"] = [date.replace('00', '') for date in split_end_date["MonthCollected2"]]
        split_end_date["DayCollected2"] = [date.replace('00', '') for date in split_end_date["DayCollected2"]]
        split_end_date["YearCollected2"] = [date.replace('0000', '') for date in split_end_date["YearCollected2"]]

        # Remove incomplete dates where day but not month was recovered
        for i in range(len(split_begin_date)):
            if split_begin_date["DayCollected"][i] != "" and split_begin_date["MonthCollected"][i] == "":
                split_begin_date["DayCollected"][i] = ""
        for i in range(len(split_end_date)):
            if split_end_date["DayCollected2"][i] != "" and split_end_date["MonthCollected2"][i] == "":
                split_end_date["DayCollected2"][i] = ""


        # Remove nonsensical dates as ambiguity
        for i in range(len(split_begin_date)):
            if split_begin_date["MonthCollected"][i] in ["02", "04", "06", "09", "11"] and split_begin_date["DayCollected"][i] == "31":
                split_begin_date["DayCollected"][i] = ""
            if split_begin_date["MonthCollected"][i] == "" and split_begin_date["DayCollected"][i] != "":
                split_begin_date["DayCollected"][i] == ""
        for i in range(len(split_end_date)):
            if split_end_date["MonthCollected2"][i] in ["02", "04", "06", "09", "11"] and split_end_date["DayCollected2"][i] == "31":
                split_end_date["DayCollected2"][i] = ""
            if split_end_date["MonthCollected2"][i] == "" and split_end_date["DayCollected2"][i] != "":
                split_end_date["DayCollected2"][i] == ""



        self.begin_date = split_begin_date
        self.end_date = split_end_date

    def prepMetadata(self):
        index = range(len(self.data))

        # Create dictionary of holding institutions
        self.essig_holdinginst['Code'] = [abbrev[0:3] for abbrev in self.essig_holdinginst['Abbrev']]
        essig_holdinginst_dict = self.essig_holdinginst.set_index('Code')['Name'].to_dict()
        
        # Create empty dataframes to hold cleaned up entries
        taxa_metadata = pd.DataFrame(index = index,\
                                      columns = ["bnhm_id", "Genus", "SpecificEpithet", "SubspecificEpithet"])
        other_metadata = pd.DataFrame(index = index,\
                                      columns = ["Taxon_Certainty", "HoldingInstitution", "EnteredBy", "BasisOfRecord", "IndividualCount"])

        # Remove jpg names
        self.data["filename"] = [fname.split("|")[0] for fname in self.data["filename"]] # new metadata in transcriptResolver function just concatenates them together, so will have to pick the first duplicate
        taxon = [name.replace('.jpg', '') for name in list(self.data["filename"])]

        # Parse out taxon uncertainty
        temp = []
        for name in taxon:
            m = re.search("_.*", name)
            if m == None:
                temp.append("")
            else:
                temp.append(m.group(0))
        
        other_metadata["Taxon_Certainty"] = temp

        # Remove underscores from taxon_certainty
        other_metadata["Taxon_Certainty"] = [taxon_cert.strip('_') for taxon_cert in other_metadata["Taxon_Certainty"]]

        # Standardize taxon certainty entries (Valid entries = ("aff.", "cf.", "gen. nr.", "near", "nr.", "poss.", "sp. nr.", "?"))
        other_metadata["Taxon_Certainty"] = [taxon_cert.replace('aff', 'aff.') for taxon_cert in other_metadata["Taxon_Certainty"]]
        other_metadata["Taxon_Certainty"] = [taxon_cert.replace('cf', 'cf.') for taxon_cert in other_metadata["Taxon_Certainty"]]
        other_metadata["Taxon_Certainty"] = [taxon_cert.replace('gen. nr', 'gen. nr.') for taxon_cert in other_metadata["Taxon_Certainty"]]
        other_metadata["Taxon_Certainty"] = [taxon_cert.replace('nr', 'nr.') for taxon_cert in other_metadata["Taxon_Certainty"]]
        other_metadata["Taxon_Certainty"] = [taxon_cert.replace('poss', 'poss.') for taxon_cert in other_metadata["Taxon_Certainty"]]
        other_metadata["Taxon_Certainty"] = [taxon_cert.replace('U', '?') for taxon_cert in other_metadata["Taxon_Certainty"]]

        # Split taxa names, but take care of possible hybrids
        taxon = [re.sub("_.*", " ", name).strip() for name in taxon]
        split_taxon = [name.split() for name in taxon]
        for x in range(len(split_taxon)):
            temp = split_taxon[x]
            if len(temp) > 4:
                if temp[3] == "x" or temp[3] == "X":
                    temp[2] = " ".join([temp[2], temp[3], temp[4]])
                    split_taxon[x] = temp[0:3]

        # Fill taxa metadata dataframe
        taxa_metadata = fill_pd(split_taxon, taxa_metadata)
        taxa_metadata['SpecificEpithet'] = [sp.replace('sp', 'sp.') for sp in taxa_metadata['SpecificEpithet']]

        # Remove filename versions (e.g., EMEC 12345.0 from bnhm_ids)
        taxa_metadata['bnhm_id'] = [re.sub("\.*", "", k) for k in taxa_metadata['bnhm_id']]

        # Dealing with subsequent changes in institution code
        taxa_metadata['bnhm_id'] = [re.sub('^CIS*', 'UCIS', k) for k in taxa_metadata['bnhm_id']]

        # Match abbreviations with respective holding institution
        other_metadata['HoldingInstitution'] = [bnhm_id[0:3] for bnhm_id in taxa_metadata["bnhm_id"]]
        other_metadata['HoldingInstitution'] = [essig_holdinginst_dict[k] for k in other_metadata['HoldingInstitution']]

        # Populate misc database fields
        other_metadata['IndividualCount'] = 1
        other_metadata['EnteredBy'] = 'Notes from Nature'
        other_metadata['BasisOfRecord'] = 'PreservedSpecimen'
        other_metadata['LifeStage'] = 'adult'
        other_metadata['PreparationType'] = 'pin'

        metadata = pd.concat([taxa_metadata, other_metadata], axis = 1) 
        self.metadata = metadata

    def normalizeGeography(self):
        index = range(len(self.data))
        
        # Create an empty dataframe to hold cleaned up entries
        split_location = pd.DataFrame(index=index, columns=["Country", "StateProvince", "County", "ContinentOcean"])

        # Ad-hoc spelling changes due to differences in Zoouniverse country list and Essig's country name standards
        self.data['Country'] = [cty.replace('Afganistan','Afghanistan') for cty in self.data['Country']]
        self.data['Country'] = [cty.replace('United Arab Erimates','United Arab Emirates') for cty in self.data['Country']]
        self.data['Country'] = [cty.replace('Iran','Iran, Islamic Republic of') for cty in self.data['Country']]
        self.data['Country'] = [cty.replace("Cote DIvoire","Cote D'Ivoire") for cty in self.data['Country']]
        self.data['Country'] = [cty.replace("Curaco","Curacao") for cty in self.data['Country']]
        self.data['Country'] = [cty.replace("Korea Sout","Korea, Republic of") for cty in self.data['Country']]
        self.data['Country'] = [cty.replace("Korea South","Korea, Republic of") for cty in self.data['Country']]
        self.data['Country'] = [cty.replace("Korea North","Korea, Democratic People's Republic of") for cty in self.data['Country']]
        self.data['Country'] = [cty.replace("Nambia","Namibia") for cty in self.data['Country']]
        self.data['Country'] = [cty.replace("Philipines","Philippines") for cty in self.data['Country']]
        self.data['Country'] = [cty.replace("Uraguay","Uruguay") for cty in self.data['Country']]
        self.data['Country'] = [cty.replace('Iran','Iran, Islamic Republic of') for cty in self.data['Country']]
        self.data['Country'] = [cty.replace('Great Britain','United Kingdom') for cty in self.data['Country']]
        self.data['Country'] = [cty.replace('Trinidad & Tobago','Trinidad and Tobago') for cty in self.data['Country']]
        self.data['Country'] = [cty.replace('Vietnam','Viet Nam') for cty in self.data['Country']]
        self.data['Country'] = [cty.replace('Antigua & Barbuda','Antigua and Barbuda') for cty in self.data['Country']]
        self.data['Country'] = [re.sub('^US$|^U\\.S\\.$|^U\\.S\\.A\\.$|^USA$|United States of America|united states','United States', cty, flags = re.IGNORECASE) for cty in self.data['Country']]

        # Normalizing country name
        # Check if country is in reference database, if there isn't an exact match, assign an NA
        print("\nChecking if country entries are valid ...")
        country_ref = list(self.essig_country['name'])
        country_index = [country_ref.index(cty)
                         if cty in country_ref else "NA"
                         for cty in self.data['Country']]

        # Populate continent ocean based on country
        print("\nPopulating continent field ...")
        continent_ref = list(self.essig_country['continent'])
        split_location["Country"] = [country_ref[index]
                                     if index is not "NA" else "NA"
                                     for index in country_index]
        split_location["ContinentOcean"] = [continent_ref[index]
                                            if index is not "NA" else "NA"
                                            for index in country_index]

        # Normalize states/province
        print("\nChecking if state and province entries are valid ...")
        state_abb_ref = list(self.essig_statecounty["State.Abbrev"])
        state_name_ref = list(self.essig_statecounty["State"])

        state_index = [
                    state_abb_ref.index(state) if state in state_abb_ref
                    else state_name_ref.index(state) if state in state_name_ref
                    else "NA"
                    for state, cty in zip(self.data['State/Province'], split_location['Country']) if cty == "United States"]

        split_location.loc[split_location['Country'] == "United States", 'StateProvince'] = [state_name_ref[index]
                                                                                    if index is not "NA" else "NA" for index in state_index]

        canada_ref = list(self.essig_canprov['name'])
        canada_index = [
                    canada_ref.index(province)
                    if province in canada_ref
                    else "NA"
                    for province, cty in zip(self.data['State/Province'], split_location['Country']) if cty == "Canada"]

        split_location.loc[split_location['Country'] == "Canada", 'StateProvince'] = [canada_ref[index]
                                                                              if index is not "NA" else "NA" for index in canada_index]


        mexico_ref = list(self.essig_mexstate['name'])
        mexico_index = [
                    mexico_ref.index(mexstate)
                    if mexstate in mexico_ref
                    else "NA"
                    for mexstate, cty in zip(self.data['State/Province'], split_location['Country']) if cty == "Mexico"]

        split_location.loc[split_location['Country'] == "Mexico", 'StateProvince'] = [mexico_ref[index]
                                                                              if index is not "NA" else "NA" for index in mexico_index]           
        # Normalize county
        print("\nChecking if county fields (only for the U.S.) are valid ...")
        county_ref = list(self.essig_statecounty['County'])

        county_index = [
            
                # Find the county in the state list
                county_ref.index(county) if county in county_ref #list(essig_statecounty.loc[essig_statecounty["State"] == state, 'County'])
                # essig_statecounty.loc[essig_statecounty["State"] == state, 'County'].index[list(essig_statecounty.loc[essig_statecounty["State"] == state, 'County']).index(county)]        
                else "NA"

                for county, state, cty in zip(self.data['County'], split_location['StateProvince'], split_location['Country']) if cty == "United States"]


        split_location.loc[split_location['Country'] == "United States", 'County'] = [county_ref[index]
                                                                                  if index is not "NA" else "NA" for index in county_index]
        # Append 'county', 'parish' and 'borough' to US counties
        split_location = split_location.fillna("")
        split_location["County"] = ["" if county == "NA" or county == ""
                                    else county + " Parish" if state == "Louisiana"
                                    else county + " Borough" if state == "Alaska"
                                    else county + " County"
                                    for county , state in zip(split_location["County"], split_location["StateProvince"])]
        
        self.geography = split_location

def main():
    ## PREAMBLE ========================
    print("\n\n\n")
    print("=" * 50)
    print("PREPARING TRANSCRIPTS FOR RESOLUTION!!")
    print("Written by Jun Ying Lim (junyinglim@gmail.com) \nfor the Essig Museum of Entomology at UC Berkeley")
    print("=" * 50)
    print("\n\n\n")
    
    ## CREATE TRANSCRIPTCLEANER INSTANCE ========================
    args = parser.parse_args()
    pipeline = transcriptCleaner(args)
    
    ## CLEANUP ========================
    ## Normalize collectors
    print("\nNormalizing collector names in the resolved transcriptions ...")
    pipeline.normalizeCollector()
    collector = pipeline.clean_collector
    
    ## Normalize dates
    print("\nNormalizing dates in the resolved transcriptions ...")
    pipeline.normalizeDates()
    
    begin_date = pipeline.begin_date
    end_date = pipeline.end_date
    
    ## Prepare metadata
    print("\nPreparing metadata fields in the resolved transcriptions ...")
    pipeline.prepMetadata()
    metadata = pipeline.metadata
    
    ## Normalize geography fields
    print("\nNormalizing geography fields in the resolved transcriptions ...")
    pipeline.normalizeGeography()
    geography = pipeline.geography

    geography["Country"]        = [x.replace("NA", "") for x in geography["Country"]]
    geography["StateProvince"]  = [x.replace("NA", "") for x in geography["StateProvince"]]
    geography["County"]         = [x.replace("NA", "") for x in geography["County"]]
    geography["ContinentOcean"] = [x.replace("NA", "") for x in geography["ContinentOcean"]]
        
    ## MERGING RESULTS ========================
    print("\nCleaning complete; merging results")
    cleanList = [metadata, collector, begin_date, end_date, geography, pipeline.data[["subject_id", "id", "Locality"]]]
    allClean = pd.concat(cleanList, axis = 1)
    allClean.rename(columns={'subject_id':'TranscriptionSubjectIDs', "id":"TranscriptionIDs"}, inplace=True)
    allClean = allClean.fillna("") # remove all NAs
    

    ## EXPORTING RESULTS ========================
    print("\nExporting results to", args.wd, "...")
    
    if args.output:
        outputfile = args.output
    else:
        outputfile = "clean_transcript"

    allClean.to_csv(os.path.join(args.wd, outputfile + ".csv"), index = False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="transcriptClean - Let's clean some transcripts!")
    parser.add_argument("-wd", help = "Working directory")
    parser.add_argument("-file", "-f", help = "File with transcriptions")
    parser.add_argument("-output", "-o", help = "Output file name")
    parser.add_argument("-username", help = "Username. Access to essig SQL database")
    parser.add_argument("-password", help = "Password. Access to essig SQL database")
    main()