from io import StringIO
from Bio import SeqIO
import Bio.SwissProt as sp
from altair.vegalite.v3.schema.channels import Key
import requests
import pandas as pd
import streamlit as st

def swissprot(uniprotID):

    url = f'https://www.uniprot.org/uniprot/{uniprotID}.txt'
    r = requests.get(url, allow_redirects=False)
    embl_str=r.text

    record_dict=SeqIO.to_dict(SeqIO.parse(StringIO(embl_str), 'swiss'))
    swiss = record_dict[list(record_dict.keys())[0]]

    #this is all wonky
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

    #pd.DataFrame([ele.split(":",1) for ele in details['comment'].split("\n")]).set_index(0).T
    comment = [ele.split(":",1) for ele in details['comment'].split("\n")]
    comment = {ele[0]:ele[1].strip() for ele in comment}
    function = comment['FUNCTION']
    try:
        biotech = comment['BIOTECHNOLOGY']
        if len(biotech) > 200:
            biotech = ""
    except KeyError:
        biotech = ""

    anno = f"{function} {biotech}"

    return pd.Series([details['name'],anno])

# def addgene(feature):
#     #Taken from bokeh_plot.py
#     #df=df.merge(featDesc,left_on="Feature",right_on="file",how="left")
#     featDesc=pd.read_csv("./feature_notes.csv",sep="\t",index_col=0)
#     hits=hits.join(featDesc)


def details(inDf):
    uniprot = inDf[inDf['uniprot']!='None']
    if not uniprot.empty:
        uniprot[["Feature","Description"]] = uniprot['uniprot'].apply(swissprot)
        uniprot['Type'] = "CDS"

    normal = inDf[inDf['uniprot']=='None']
    featDesc=pd.read_csv("./feature_notes.csv",sep="\t")
    normal=normal.merge(featDesc,left_on = "sseqid", right_on = "file", how = "left")
    
    outDf = uniprot.append(normal)

    return outDf
