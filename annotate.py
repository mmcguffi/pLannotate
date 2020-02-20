import numpy as np
import base64

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from tempfile import NamedTemporaryFile
import pandas as pd
import streamlit as st

def BLAST(seq,wordsize=12, db='nr_db', BLASTtype="p", flags = 'qstart qend sseqid sframe pident slen sseq length sstart send qlen'):
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    subprocess.call(
        (f'blast{BLASTtype} -task blastn-short -query {query.name} -out {tmp.name} '
        f'-db {db} -max_target_seqs 20000 -word_size {str(wordsize)} -outfmt "6 {flags}"'),
        shell=True)

    with open(tmp.name, "r") as file_handle:  #opens BLAST file
        align = file_handle.readlines()

    tmp.close()
    query.close()

    inDf=pd.DataFrame([ele.split() for ele in align],columns=flags.split())
    inDf=inDf.apply(pd.to_numeric, errors='ignore')

    return inDf
def clean_and_calculate(inDf):
    inDf['qstart']=inDf['qstart'].astype('int')
    inDf['qstart']=inDf['qstart']-1
    inDf['qend']=inDf['qend']-1
    inDf['percmatch'] = (inDf['length']/inDf['slen']*100)
    inDf[['sseqid','type']]=inDf['sseqid'].str.split("|", n=1, expand=True)
    inDf['sseqid']=inDf['sseqid'].str.replace(".gb","")#gb artifact from BLASTdb
    inDf['abs percmatch']=100-abs(100-inDf['percmatch'])#eg changes 102.1->97.9
    inDf['pi_permatch']=(inDf["pident"]*inDf["abs percmatch"])/100
    inDf['score']=(inDf['pi_permatch']/100)*inDf["length"]
    inDf['qlen']=(inDf['qlen']/2).astype('int')
    inDf['fragment'] = inDf["percmatch"] < 95
    inDf['Feature'] = inDf['sseqid'].str.replace("_"," ") #idk where this happened in other scripts

    inDf=inDf.sort_values(by=["score","length","percmatch"], ascending=[False, False, False])

    #subtracts a full plasLen from starts > plas_len
    #subtracts a full plasLen from ends > plas_len IF start ALSO > plasLen
    inDf['qend']=np.where((inDf['qend']>=inDf['qlen']) & (inDf['qstart']>=inDf['qlen']),inDf['qend']-inDf['qlen'],inDf['qend'])
    inDf['qstart']=np.where(inDf['qstart']>=inDf['qlen'],inDf['qstart']-inDf['qlen'],inDf['qstart'])
    inDf=inDf.drop_duplicates()

    #remove overlapping hits -- need to add wiggle, fix origin
    #goes from bottom of list and looks if it lies within
    #the any of the fragements above, which are sorted best
    #to worst
    dropIndexes=[]
    inDf=inDf.reset_index(drop=True)
    wiggle=6 #this also can through lists out of index in the current way it works
             #maybe change to a percentage of overlap instead of absolute bps?
    inDf['wstart']=inDf['qstart']+wiggle
    inDf['wend']=inDf['qend']+(wiggle-1)

    #this loop is slow -- make it a single pandas filter
    #can acheive this by `index > curIndex`? something like that
    for i in range(len(inDf)-1,0,-1):
        start=(inDf.iloc[i]['wstart'] >= inDf.iloc[list(range(0,i))][['wstart']]).wstart & (inDf.iloc[i]['wstart'] <= inDf.iloc[list(range(0,i))][['wend']]).wend
        end=(inDf.iloc[i]['wend'] >= inDf.iloc[list(range(0,i))][['wstart']]).wstart & (inDf.iloc[i]['wend'] <= inDf.iloc[list(range(0,i))][['wend']]).wend
        if True in (start|end).values: dropIndexes.append(i)
    inDf=inDf.drop(dropIndexes)

    #subtracts a full plasLen from ends > plas_len
    inDf['qend']=np.where(inDf['qend']>=inDf['qlen'] ,inDf['qend']-inDf['qlen'],inDf['qend'])
    inDf=inDf.astype({'qstart': 'int','qend': 'int'})

    return inDf
def FeatureLocation_smart(r):
    #creates compound locations if needed
    if r.qend>r.qstart:
        return FeatureLocation(r.qstart, r.qend, r.sframe)
    elif r.qstart>r.qend:
        first=FeatureLocation(r.qstart, r.qlen, r.sframe)
        second=FeatureLocation(0, r.qend, r.sframe)
        if r.sframe == 1 or r.sframe == 0:
            return first+second
        elif r.sframe == -1:
            return second+first
def get_gbk(inDf,inSeq):
    #adds a FeatureLocation object so it can be used in gbk construction
    inDf['feat loc']=inDf.apply(FeatureLocation_smart, axis=1)

    #make a record
    record = SeqRecord(seq=Seq(inSeq),name='pLannotate')
    record.seq.alphabet=generic_dna
    record.annotations["topology"] = "circular"
    for index in inDf.index:
        record.features.append(SeqFeature(inDf.loc[index]['feat loc'], type=inDf.loc[index]["type"],qualifiers={"label": index, "identity":inDf.loc[index]["pident"],"match length":inDf.loc[index]["percmatch"], "Other:":inDf.loc[index]["type"]}))

    #converts gbk into straight text
    outfileloc=NamedTemporaryFile()
    with open(outfileloc.name, "w") as handle:
        SeqIO.write(record, handle, "genbank")
    with open(outfileloc.name) as handle:
        record=handle.read()
    outfileloc.close()

    return record

def annotate(inSeq):
    database="./BLAST_dbs/full_snapgene_feature_list_w_types_db"

    #I could just create a seq object? this could catch errors though
    fileloc = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(inSeq),name="pLannotate"), fileloc.name, 'fasta')
    record=list(SeqIO.parse(fileloc.name, "fasta"))
    fileloc.close()

    assert len(record)==1,f"FASTA file contains ~multitudes~ --> please submit a single FASTA file."
    record=record[0]

    query=str(record.seq)*2

    blastDf = BLAST(seq=query,wordsize=12, db=database, BLASTtype="n")
    if blastDf.empty: #if no hits are found
        return blastDf
    blastDf = clean_and_calculate(blastDf)

    smallHits=blastDf[blastDf['slen']<25]
    smallHits=smallHits[smallHits["pident"] >= ((smallHits["slen"]-1)/smallHits["slen"])*100] #allows for 1 mismatch
    smallHits=smallHits[smallHits["percmatch"] >= ((smallHits["slen"]-1)/smallHits["slen"])*100]
    smallHits['small']=True

    normHits=blastDf[blastDf['slen']>=25]
    normHits=normHits[normHits['length']>=18]
    normHits=normHits[normHits['pident']>=95]
    normHits['small']=False

    hits=smallHits.append(normHits)
    hits=hits.set_index(['Feature'])

    st.write("all hits")
    st.write(hits)

    return hits
