import subprocess
from tempfile import NamedTemporaryFile

import pandas as pd
import streamlit as st

import plannotate.resources as rsc

def details(inDf, blast_database):
    
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
    databases = rsc.get_yaml(blast_database)
    
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
            
            details_file_loc = db_details['location']
            #if the description file is compressed
            if rsc.get_name_ext(details_file_loc)[1] == '.gz':
                feat_desc = parse_gz(sseqids, details_file_loc)
            else: #if it is uncompressed
                feat_desc = pd.read_csv(details_file_loc)
                #featDesc['sseqid'] = featDesc['sseqid'].astype(str)
        
        #if no file is passed, data should already be in dataframe
        else:
            feat_desc = inDf.loc[inDf['db'] == database_name][['sseqid','Feature','Description']]
            
        #try to see if a default type was passed
        try:
            default_type = db_details['default_type']
            feat_desc['Type'] = default_type
        except KeyError:
            pass 
        
        details_list.append(feat_desc)
        
    details_list = pd.concat(details_list)
    
    #try dropping extra 'Feature' and 'Description' columns so it's not duplicated
    #if it already existed it should be saved above
    try:
        inDf = inDf.drop("Description", axis = 1)
        inDf = inDf.drop("Feature", axis = 1)
    except KeyError:
        pass
    
    out_df = inDf.merge(details_list, on='sseqid', how='left')
        

    # uniprot = inDf[inDf['db']=='swissprot'].copy()
    # if not uniprot.empty:

    #     hits = uniprot['sseqid'].tolist()
    #     hits = [_ for _ in hits if _] #removes blank edgecases
    #     hits = "|".join(hits)
    #     output = NamedTemporaryFile(suffix="csv")
    #     subprocess.call(f'rg -z "{hits}" {rsc.get_resource("data","swiss_description_verbose.csv.gz")} > {output.name}',shell = True)
    #     df = pd.read_csv(output.name, header = None, names=['sseqid','Feature','Description'])
    #     df['Type'] = "swissprot"
    #     output.close()
    #     uniprot = uniprot.merge(df, how = 'left', on = 'sseqid')

    # fpbase = inDf[inDf['db']=='fpbase'].copy()
    # if not fpbase.empty:
    #     fpblurb = pd.read_csv(rsc.get_resource("data", "fpbase_burbs.csv"),index_col=0)
    #     fpbase['Type'] = "CDS" # uppercase is for coloring (colors.csv)
    #     fpbase['Feature'] = fpbase['sseqid']
    #     fpbase = fpbase.merge(fpblurb, on = "sseqid", how = 'left')

    # infernal = inDf[inDf['db']=='infernal'].copy()
    # if not infernal.empty:
    #     infernal['Type'] = "ncRNA" # uppercase is for coloring (colors.csv)
    #     # infernal = infernal.rename(columns = {"target name":"Feature","description of target":"Description"})
    #     # infernal['Feature'] = infernal['Feature'].str.replace("_"," ")
    #     # infernal['Description'] = "Accession: " + infernal['accession'] + " - " + infernal['Description']

    # addgene = inDf[inDf['db']=='addgene'].copy()
    # featDesc=pd.read_csv(rsc.get_resource("data", "addgene_collected_features_test_20-12-11_description.csv"))
    # addgene=addgene.merge(featDesc, on = "sseqid", how = "left")

    # outDf = uniprot.append(addgene)
    # outDf = outDf.append(fpbase)
    # outDf = outDf.append(infernal)

    return out_df
