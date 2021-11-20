import subprocess
from tempfile import NamedTemporaryFile

import pandas as pd
import streamlit as st

import plannotate.resources as rsc

def details(inDf):
    
    def parse_gz(sseqids, gz_loc):
    #this is a bit fragile right now -- requires ['sseqid','Feature','Description'] order
    #as well as a default type
    #currently this is only implemented for the large SwissProt db
    #Could scrape first line to infer what is given that way
        hits = "|".join(sseqids)
        output = NamedTemporaryFile(suffix="csv")
        subprocess.call(f'rg -z "{hits}" {gz_loc} > {output.name}',shell = True)
        gz_details = pd.read_csv(output.name, header = None, names=['sseqid','Feature','Description'])
        output.close()
        return gz_details
    
    #loop through databases
    databases = rsc.get_yaml()
    
    dbs_used = set(inDf['db'].to_list())
    
    details_list = []
    for database_name in dbs_used:
        database = databases[database_name]
        
        sseqids = inDf.loc[inDf['db'] == database_name]['sseqid'].tolist()
        sseqids = [_ for _ in sseqids if _] #removes blank edgecases

        db_details = database['details']
        if database['method'] == 'infernal':
            pass

        if db_details['file'] == True:
            
            if db_details['location'] == "Default":
                details_file_loc = rsc.get_details(database_name) + ".csv"
            else:
                details_file_loc = db_details['location']

            #if the description file is compressed
            if db_details['compressed'] == True:
                if rsc.get_name_ext(details_file_loc)[1] != ".gz":
                    details_file_loc +=  ".gz"
                feat_desc = parse_gz(sseqids, details_file_loc) 
                
            else: #if it is uncompressed
                feat_desc = pd.read_csv(details_file_loc)
            
            # bespoke extraction of swissprot protein exisitence level
            if database_name == 'swissprot':
                level = feat_desc['Description'].str.find("existence level") + 16 #len of "existence level" + 1
                feat_desc['s'] = level
                feat_desc['e'] = level + 1

                def calc_priority_mod(d,s,e):
                    if s == 15 and e == 16:
                        return 0
                    else:
                        return int(d[s:e]) - 1
    
                # extract the level from the description
                feat_desc['priority_mod'] = [calc_priority_mod(d,s,e) for d, s, e in zip(feat_desc["Description"], feat_desc["s"], feat_desc["e"])]

                feat_desc = feat_desc.drop(columns=['s','e'])

        #if no file is passed, data should already be in dataframe
        else:
            feat_desc = inDf.loc[inDf['db'] == database_name][['sseqid','Feature','Description']]
            
        #try to see if a default type was passed
        try:
            default_type = db_details['default_type']
            feat_desc['Type'] = default_type
        except KeyError:
            pass 
        
        #details_list.append(feat_desc)
        
    #details_list = pd.concat(details_list)
    
    #try dropping extra 'Feature' and 'Description' columns so it's not duplicated
    #if it already existed it should be saved above
    try:
        inDf = inDf.drop("Description", axis = 1)
        inDf = inDf.drop("Feature", axis = 1)
    except KeyError:
        pass
    
    #if not details_list.empty:
    out_df = inDf.merge(feat_desc, on='sseqid', how='left')
    #else:
    #    out_df = pd.DataFrame()
    
    #drop primers -- not super useful
    # st.write(out_df)   
    # out_df = out_df[out_df['Type'] != 'primer bind']
    # st.write(out_df)   
    return out_df
