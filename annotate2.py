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

def FeatureLocation_smart(r):
    if r.qend>r.qstart:
        return FeatureLocation(r.qstart, r.qend, r.sframe)
    elif r.qstart>r.qend:
        first=FeatureLocation(r.qstart, r.qlen, r.sframe)
        second=FeatureLocation(0, r.qend, r.sframe)
        if r.sframe == 1 or r.sframe == 0:
            return first+second
        elif r.sframe == -1:
            return second+first

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

    alignDf=pd.DataFrame([ele.split() for ele in align],columns=flags.split())
    alignDf=alignDf.apply(pd.to_numeric, errors='ignore')
    alignDf['qstart']=alignDf['qstart'].astype('int')
    alignDf['qstart']=alignDf['qstart']-1
    alignDf['qend']=alignDf['qend']-1
    alignDf['percmatch'] = (alignDf['length']/alignDf['slen']*100)
    alignDf[['sseqid','type']]=alignDf['sseqid'].str.split("|", n=1, expand=True)
    alignDf['sseqid']=alignDf['sseqid'].str.replace(".gb","")
    alignDf['abs percmatch']=100-abs(100-alignDf['percmatch'])#eg changes 102.1->97.9
    alignDf['pi_permatch']=(alignDf["pident"]*alignDf["abs percmatch"])/100
    alignDf['score']=(alignDf['pi_permatch']/100)*alignDf["length"]
    alignDf['qlen']=(alignDf['qlen']/2).astype('int')
    #alignDf['plas_len']=len(seq)

    alignDf=alignDf.sort_values(by=["score","length","percmatch"], ascending=[False, False, False])
    # alignDf=alignDf.drop(alignDf[(alignDf['qstart']>=alignDf['qlen']/2)&(alignDf['qend']>=alignDf['qlen']/2)].index)

    alignDf['qstart']=np.where(alignDf['qstart']>=len(seq)/2,alignDf['qstart']-len(seq)/2,alignDf['qstart'])
    alignDf['qend']=np.where(alignDf['qend']>=len(seq)/2,alignDf['qend']-len(seq)/2,alignDf['qend'])
    alignDf=alignDf.drop_duplicates()
    alignDf=alignDf.astype({'qstart': 'int','qend': 'int'})

    #alignDf.to_csv("~/Desktop/test.csv")

    alignDf['feat loc']=alignDf.apply(FeatureLocation_smart, axis=1)
    st.write(alignDf)

    return align, alignDf

def get_hits(inHits):
    df = []

    for ele in inHits:
        splitAlign=ele.split()

        qstart = int(splitAlign[0])-1
        qend = int(splitAlign[1])-1
        sframe = int(splitAlign[3])
        pident = float(splitAlign[4])
        slen = int(splitAlign[5])
        sseq = splitAlign[6]
        percmatch = round((len(sseq)/slen)*100,3)

        sseqid = splitAlign[2]
        partType=sseqid.split("|")[1]
        sseqid=sseqid.split("|")[0]
        name = sseqid.split(".gb")[0].replace("_"," ")#.split(" (")[0]

        abspercmatch=100-abs(100-percmatch)#eg changes 102.1->97.9

        absdiff=(pident/100)*(abspercmatch/100)*len(sseq)
        df.append({'Abs. diff': absdiff, 'name': name,'type':partType,'start':qstart,'end':qend,'frame':sframe, 'percent identity': pident, 'percent match': abspercmatch, "Length of hit":len(sseq),"Length of target seq":slen} )
    df=pd.DataFrame(df)

    st.write("current one")
    st.write(df)#########

    df=df.sort_values(by=["Abs. diff","Length of hit",'percent match'], ascending=[False, False, False])
    chunk=df[df["type"]=='source']
    #df=df[df["type"]!='source']

    return df,chunk

def annotate(inSeq):
    database="./BLAST_dbs/full_snapgene_feature_list_w_types_db"

    recordDf=pd.DataFrame()
    wiggle=6 #maybe change to a percentage of overlap instead of absolute bps?
            #this also can through lists out of index in the current way it works

    fileloc = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(inSeq),name="Temp"), fileloc.name, 'fasta')

    record=list(SeqIO.parse(fileloc.name, "fasta"))
    fileloc.close()
    assert len(record)==1,f"FASTA file contains ~multitudes~ --> please submit a single FASTA file."
    record=record[0]
    record.seq.alphabet=generic_dna
    record.annotations["topology"] = "circular"
    query=str(record.seq)*2
    pLen = len(str(record.seq))
    #st.write(pLen)
    seqSpace=[[] for i in range(len(query))]
    #first annotates full small features, 12-25 nts
    blast,blastDf =BLAST(seq=query,wordsize=12, db=database, BLASTtype="n")

    #if blast.empty:#for df
    if not blast:
        return None, None

    else:
        hits,chunk = get_hits(blast)
        #hits=blast

        smallHits=blastDf[blastDf['slen']<25]
        smallHits=smallHits[smallHits["pident"] >= ((smallHits["slen"]-1)/smallHits["slen"])*100] #allows for 1 mismatch
        smallHits=smallHits[smallHits["percmatch"] >= ((smallHits["slen"]-1)/smallHits["slen"])*100]

        normHits=blastDf[blastDf['slen']>=25]

        st.write("small hits")
        st.write(smallHits)

        for ele in hits.index:
            slen=int(hits.loc[[ele]]['Length of target seq'])

            if slen < 25:
                qstart = int(hits.loc[[ele]]['start'])
                qend = int(hits.loc[[ele]]['end'])
                sframe = int(hits.loc[[ele]]['frame'])
                pident = float(hits.loc[[ele]]['percent identity'])
                percmatch = float(hits.loc[[ele]]['percent match'])
                name = hits.loc[[ele]]['name'].values[0]
                partType=hits.loc[[ele]]['type'].values[0]

                #if pident==100 and percmatch==100: #filters out all other hits
                if pident >= ((slen-1)/slen)*100: #filters out all other hits
                    if percmatch >= ((slen-1)/slen)*100: #length of match

                        if qend-pLen>0 and qstart-pLen>=0:
                            continue
                        elif qend-pLen>0 and qstart-pLen<=0:
                            first=FeatureLocation(qstart, pLen,sframe)
                            second=FeatureLocation(0, qend-pLen+1,sframe)
                            if sframe == 1 or sframe == 0:
                                featLoc=first+second
                            elif sframe == -1:
                                featLoc=second+first
                            else:
                                print("error1.1")
                        elif qend<=pLen and qstart<=pLen:
                            featLoc=FeatureLocation(qstart, qend+1,sframe)
                        else:
                            print("error1.2")

                        record.features.append(SeqFeature(featLoc, type=partType,qualifiers={"label": name,"identity":pident,"match length":percmatch, "Other:":partType}))
                        recordDf=recordDf.append(hits.loc[[ele]])

                        for i in featLoc:
                            seqSpace[i].append((name,sframe,pident,percmatch))
                            if len(featLoc.parts) > 1:
                                seqSpace[i+len(record.seq)].append((name,sframe,pident,percmatch)) #unique id -- necessary?

        blast,blastDf=BLAST(seq=query,wordsize=18, db =database,BLASTtype="n")
        hits,chunk = get_hits(blast)

        #st.write(hits)
        for ele in hits.index:
            qstart = int(hits.loc[[ele]]['start'])
            qend = int(hits.loc[[ele]]['end'])

            seqSpaceSlice=seqSpace[qstart+wiggle:qend-(wiggle-1)]

            seqSpaceSlice=[set(ele) for ele in seqSpaceSlice]
            occupiedSpace=set.intersection(*seqSpaceSlice)

            if not occupiedSpace:
                name = hits.loc[[ele]]['name'].values[0]
                sframe = int(hits.loc[[ele]]['frame'])
                pident = float(hits.loc[[ele]]['percent identity'])
                percmatch = float(hits.loc[[ele]]['percent match'])
                partType=hits.loc[[ele]]['type'].values[0]

                if pident>=95: #filters out all other hits

                    if qend-pLen>0 and qstart-pLen>=0:
                        continue
                    elif qend-pLen>0 and qstart-pLen<=0:
                        first=FeatureLocation(qstart,pLen,sframe)
                        second=FeatureLocation(0, qend-pLen+1,sframe)
                        if sframe == 1 or sframe == 0:
                            featLoc=first+second
                        elif sframe == -1:
                            featLoc=second+first
                        else:
                            print("error2.1")
                    elif qend<=pLen and qstart<=pLen:
                        featLoc=FeatureLocation(qstart, qend+1,sframe)
                    else:
                        print("error2.2")

                    if percmatch < 95: #length of match
                        Type="Fragment"
                    else:
                        Type=str(partType)

                    for i in featLoc:
                        seqSpace[i].append((name,sframe,pident,percmatch))
                        if len(featLoc.parts) > 1:
                            seqSpace[i+len(record.seq)].append((name,sframe,pident,percmatch))

                    #if name != "ColE1 ori chunk":
                    record.features.append(SeqFeature(featLoc, type=Type,qualifiers={"label": name,"identity":pident,"match length":percmatch, "Other:":partType}))
                    recordDf=recordDf.append(hits.loc[[ele]])

        record.name=record.name[:20].replace(" ","_")#######maybe change this

        # if fragmentMode == False:
        # 	wholeFeats=[]
        # 	for feat in record.features:
        # 		if feat.type != "Fragment" or feat.qualifiers['label'] == "ColE1 ori chunk":
        # 			wholeFeats.append(feat)
        # 	record.features=wholeFeats

        recordDf=recordDf.sort_values(by=["Abs. diff"],ascending=[False])
        recordDf=recordDf.drop("Abs. diff",axis=1).set_index("name",drop=True)


        return record, recordDf
