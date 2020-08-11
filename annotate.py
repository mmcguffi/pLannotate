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
        (f'blast{BLASTtype} -task blastn-short -query {query.name} -out {tmp.name} -perc_identity 95 ' #pi needed?
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

    #Activity selection problem -- greedy -- possible solution?
    #https://iq.opengenus.org/activity-selection-problem-greedy-algorithm/
    #interval tree package
    #pandas.Interval.overlaps
    #pandas.IntervalIndex.contains

    #remove overlapping hits -- need to add wiggle, fix origin
    #goes from bottom of list and looks if it lies within
    #the any of the fragements above, which are sorted best
    #to worst
    dropIndexes=[]
    inDf=inDf.reset_index(drop=True)
    wiggle=7 #this also can through lists out of index in the current way it works
             #maybe change to a percentage of overlap instead of absolute bps?
    inDf['wstart']=inDf['qstart']+wiggle
    inDf['wend']=inDf['qend']-(wiggle-1)

    #this loop is slow -- make it a single pandas filter
    #can acheive this by `index > curIndex`? something like that
    st.write("raw hits")
    st.write(inDf)
    import time
    stT=time.time()
    #review code below
    inDf['qend']=np.where(inDf['qstart']==0, inDf['qend']-inDf['qlen'], inDf['qend'])
    inDf['qstart']=np.where(inDf['qstart']==0, inDf['qstart']-inDf['qlen'], inDf['qstart'])

    # if one is > plas len, subtract plaslen
    # inDf['qlen_copy']=inDf['qlen']
    # my_cols = ['qstart','qend','qlen_copy']
    # inDf[['qstart','qend']]=inDf[my_cols].where(~(inDf['qend']>inDf['qlen_copy'])|(inDf['qstart']>inDf['qlen_copy']),inDf[my_cols].apply(lambda x: x-x['qlen_copy'],axis=1)).drop("qlen_copy",axis=1)
    # inDf

    #inDf[['qstart','qend']]=np.where((inDf['qend']>inDf['qlen'])|(inDf['qstart']>inDf['qlen']), inDf[inDf['qstart']-inDf['qlen'],inDf['qend']-inDf['qlen']], inDf[['qstart','qend']])
    inDf['qstart']=np.where(inDf['qstart']>inDf['qlen'], inDf['qstart']-inDf['qlen'], inDf['qstart'])

    # #################################################################
    # for i in range(len(inDf) - 1, 0, -1):
    #     start=(inDf.iloc[i]['wstart'] >= inDf.iloc[list(range(0,i))][['wstart']]).wstart & (
    #            inDf.iloc[i]['wstart'] <= inDf.iloc[list(range(0,i))][['wend']]).wend
    #
    #     end  =(inDf.iloc[i]['wend']   >= inDf.iloc[list(range(0,i))][['wstart']]).wstart & (
    #            inDf.iloc[i]['wend']   <= inDf.iloc[list(range(0,i))][['wend']]).wend
    #
    #     if True in (start|end).values: dropIndexes.append(i)
    # inDf=inDf.drop(dropIndexes)
    # #################################################################

    l=[]
    inDf['drop']=False
    inDf['level']=0
    inDf=inDf[inDf['pident']>=95]
    for i in inDf.index:
        df=inDf[inDf.index<i]
        df=df[df['drop']==False]
        s=inDf.loc[i]['wstart']
        e=inDf.loc[i]['wend']
        # #####
        # # if i == 28:
        # st.write(inDf.loc[i]['sseqid'],i,df)
        # st.write(s,e)
        # srt=((df['wstart']<=s) & (df['wend']>=s))
        # en=((df['wstart']<=e) & (df['wend']>=e))
        # st.write(df[srt&en])
        # st.write("===")
        # #####
        startBound=((df['qstart']<=s) & (df['qend']>=s))
        endBound=((df['qstart']<=e) & (df['qend']>=e))
        if df[startBound].empty ^ df[endBound].empty: level=1 # ^ == XOR
        else: level=0
        within=df[startBound&endBound]
        #st.write(within)
        if not within.empty: drop=True
        else: drop=False
        #l.append((level,drop))
        inDf.at[i,'drop'] = drop
        inDf.at[i,'level'] = level
    #l=pd.DataFrame(l,columns=['level','drop'])
    inDf=inDf[inDf['drop']==False]
    st.write("dropped",inDf)

    st.write(time.time()-stT)

    #subtracts a full plasLen from ends > plas_len
    inDf['qend']=np.where(inDf['qend']>=inDf['qlen'] ,inDf['qend']-inDf['qlen'],inDf['qend'])
    inDf['qstart']=np.where(inDf['qstart']>=inDf['qlen'] ,inDf['qstart']-inDf['qlen'],inDf['qstart'])
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
    #this could be passed a more annotated df
    inDf=inDf.reset_index(drop=True)

    #adds a FeatureLocation object so it can be used in gbk construction
    inDf['feat loc']=inDf.apply(FeatureLocation_smart, axis=1)
    #make a record
    record = SeqRecord(seq=Seq(inSeq),name='pLannotate')
    record.seq.alphabet=generic_dna
    record.annotations["topology"] = "circular"
    for index in inDf.index:
        record.features.append(SeqFeature(
            inDf.loc[index]['feat loc'],
            type=inDf.loc[index]["type"],
            qualifiers={"label": inDf.loc[index]["sseqid"].replace("_"," "),
            "identity":inDf.loc[index]["pident"],
            "match length":inDf.loc[index]["percmatch"],
            "Other:":inDf.loc[index]["type"]}))

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

    return hits
