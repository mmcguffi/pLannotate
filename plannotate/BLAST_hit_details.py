from io import StringIO
import re
import requests

from Bio import SeqIO
import Bio.SwissProt as sp
# Seems unused
#from altair.vegalite.v3.schema.channels import Key
import pandas as pd
import streamlit as st

import plannotate

def swissprot(uniprotID):

    url = f'https://www.uniprot.org/uniprot/{uniprotID}.txt'
    r = requests.get(url, allow_redirects=False)
    embl_str=r.text

    record_dict=SeqIO.to_dict(SeqIO.parse(StringIO(embl_str), 'swiss'))
    swiss = record_dict[list(record_dict.keys())[0]]

    #this is all wonky and probably over-engineered
    attrs = ["annotations,gene_name","description","annotations,comment","annotations,protein_existence"]
    details = {"name":swiss.name}
    for attr in attrs:
        attr=attr.split(",")
        try:
            if len(attr)==1:
                details[attr[1]] = getattr(swiss,attr[0])
            elif len(attr)==2:
                details[attr[1]] = getattr(swiss,attr[0])[attr[1]]
        except:
            pass
        
    try:    
        protein_existence = pd.read_csv("./data/protein_existence.csv",header=None,index_col = 0, names= ['description'])
        protein_existence = protein_existence.loc[details['protein_existence']]['description'] + ". "
    except (KeyError, IndexError):
        protein_existence = ""

    #pd.DataFrame([ele.split(":",1) for ele in details['comment'].split("\n")]).set_index(0).T
    try:
        comment = [ele.split(":",1) for ele in details['comment'].split("\n")]
        comment = {ele[0]:ele[1].strip() for ele in comment}
    except KeyError:
        comment = {}

    try:
        function = comment['FUNCTION']
        function = f"{function} " #gives a space for stuff after
    except KeyError:
        function = ""

    try:
        biotech = comment['BIOTECHNOLOGY']
        #cuts out super long details (eg GFP)
        if len(biotech) > 200:
            biotech = ""
    except KeyError:
        biotech = ""

    try:
        organism = swiss.annotations['organism']
        organism = f"From {organism}. "
    except KeyError:
        organism = ""

    #tries to get a common name and an alt name, if available

    try:
        names = details['gene_name'].split(";")

        name = [ele for ele in names if "Name=" in ele][0]
        name = name.replace("Name=","")

        swissprotName = swiss.name
        swissprotName = f"{swissprotName} - " #gives a space for stuff after

    except (KeyError, IndexError):
        #altName = None
        name = details['name']
        swissprotName = ""


    try:
        names = details['gene_name'].split(";")
        altName = [ele for ele in names if "Synonyms=" in ele][0]
        altName = altName.replace("Synonyms=","")
        altName = f"Also known as {altName}. "
    except (KeyError, IndexError):
        altName = ""


    anno = f"{swissprotName}{altName}{protein_existence}{function}{organism}{biotech}"
    anno = re.sub(r"\s*{.*}\s*", " ", anno).strip() #removes text between curly braces
    name = re.sub(r"\s*{.*}\s*", " ", name).strip() #this removes pubmed citations and
                                                    #extra details on some names
                                            
    return pd.Series([name, anno])

def details(inDf):

    uniprot = inDf[inDf['db']=='swissprot'].copy()
    if not uniprot.empty:
        uniprot[["Feature","Description"]] = uniprot['uniprot'].apply(swissprot)
        uniprot['Type'] = "swissprot" # for coloring (colors.csv)

    fpbase = inDf[inDf['db']=='fpbase'].copy()
    if not fpbase.empty:
        fpblurb = pd.read_csv(plannotate.get_resource("data", "fpbase_burbs.csv"),index_col=0)
        fpbase['Type'] = "CDS" # uppercase is for coloring (colors.csv)
        fpbase['Feature'] = fpbase['sseqid']
        fpbase = fpbase.merge(fpblurb, on = "sseqid", how = 'left')

    infernal = inDf[inDf['db']=='infernal'].copy()
    if not infernal.empty:
        infernal['Type'] = "ncRNA" # uppercase is for coloring (colors.csv)
        infernal = infernal.rename(columns = {"target name":"Feature","description of target":"Description"})
        infernal['Feature'] = infernal['Feature'].str.replace("_"," ")
        infernal['Description'] = "Accession: " + infernal['accession'] + " - " + infernal['Description']

    addgene = inDf[inDf['db']=='addgene'].copy()
    featDesc=pd.read_csv(plannotate.get_resource("data", "addgene_collected_features_test_20-12-11_description.csv"))
    addgene=addgene.merge(featDesc, on = "sseqid", how = "left")

    outDf = uniprot.append(addgene)
    outDf = outDf.append(fpbase)
    outDf = outDf.append(infernal)

    return outDf
