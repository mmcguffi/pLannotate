import subprocess
from tempfile import NamedTemporaryFile

import pandas as pd
import streamlit as st

import plannotate

def details(inDf):

    uniprot = inDf[inDf['db']=='swissprot'].copy()
    if not uniprot.empty:

        hits = uniprot['sseqid'].tolist()
        hits = [_ for _ in hits if _] #removes blank edgecases
        hits = "|".join(hits)
        output = NamedTemporaryFile(suffix="csv")
        subprocess.call(f'rg -z "{hits}" {plannotate.get_resource("data","swiss_description_verbose.csv.gz")} > {output.name}',shell = True)
        df = pd.read_csv(output.name, header = None, names=['sseqid','Feature','Description'])
        df['Type'] = "swissprot"
        output.close()
        uniprot = uniprot.merge(df, how = 'left', on = 'sseqid')

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
