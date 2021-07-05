import base64
import os
import subprocess
from tempfile import NamedTemporaryFile

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import numpy as np
import pandas as pd
import streamlit as st

from plannotate.infernal import parse_infernal


def BLAST(seq,wordsize=12, db='nr_db', task="BLAST"):
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    if task == "BLAST":
        flags = 'qstart qend sseqid sframe pident slen sseq length sstart send qlen evalue'
        subprocess.call( #remove -task blastn-short?
            (f'blastn -task blastn-short -query {query.name} -out {tmp.name} -perc_identity 95 ' #pi needed?
            f'-db {db} -max_target_seqs 20000 -culling_limit 25 -word_size {str(wordsize)} -outfmt "6 {flags}"'),
            shell=True)

    elif task == "DIAMOND":
        flags = 'qstart qend sseqid pident slen length sstart send qlen evalue'
        extras = '-l 1 --matrix PAM30 --id 10 --quiet'
        subprocess.call(f'diamond blastx -d {db} -q {query.name} -o {tmp.name} '
                        f'{extras} --outfmt 6 {flags}',shell=True)

    elif task == "infernal":
        flags = "--cut_ga --rfam --nohmmonly --fmt 2"
        cmd = f"cmscan {flags} --tblout {tmp.name} --clanin {db} {query.name}"
        print(cmd)
        subprocess.call(cmd, shell=True)

        inDf = parse_infernal(tmp.name)

        tmp.close()
        query.close()

        return inDf

    with open(tmp.name, "r") as file_handle:  #opens BLAST file
        align = file_handle.readlines()

    tmp.close()
    query.close()

    inDf = pd.DataFrame([ele.split() for ele in align],columns=flags.split())
    inDf = inDf.apply(pd.to_numeric, errors='ignore')

    return inDf

def calc_level(inDf):
    #calculates the level to be rendered at
    #this must happen while the concatentated seq df still exists
    inDf['level']=0
    for i in inDf.index:
        df=inDf[inDf.index<i]
        s=inDf.loc[i]['qstart']
        e=inDf.loc[i]['qend']
        startBound=((df['qstart']<=s) & (df['qend']>=s))
        endBound=((df['qstart']<=e) & (df['qend']>=e))
        #st.write(startBound.sum())
        if df[startBound].empty ^ df[endBound].empty:
            level=1 # ^ == XOR
        else:
            level=0
        within=df[startBound&endBound]
        #st.write(within)
        inDf.at[i,'level'] = level
    return inDf

def calculate(inDf, task, is_linear):

    inDf['qstart'] = inDf['qstart']-1
    inDf['qend']   = inDf['qend']-1

    if task == "BLAST":
        inDf['uniprot'] = 'None'
        inDf['priority'] = 0

    elif task == "DIAMOND":
        try:
            inDf[['sp','uniprot','sseqid']] = inDf['sseqid'].str.split("|", n=2, expand=True)
        except ValueError:
            pass
        inDf['sframe'] = (inDf['qstart']<inDf['qend']).astype(int).replace(0,-1)
        inDf['slen']   = inDf['slen'] * 3
        inDf['length'] = abs(inDf['qend']-inDf['qstart'])+1
        inDf['priority'] = 1

    elif task == "infernal":
        inDf["priority"] = 2
        inDf['uniprot'] = 'None'
        inDf['sseq'] = ""
        inDf["sframe"] = inDf["sframe"].replace(["-","+"], [-1,1])
        inDf['qstart'] = inDf['qstart']-1
        inDf['qend']   = inDf['qend']-1
        inDf['length'] = abs(inDf['qend']-inDf['qstart'])+1
        inDf['slen'] = abs(inDf['send']-inDf['sstart'])+1
        inDf['pident'] = 100

    inDf = inDf[inDf['evalue'] < 1].copy() #gets rid of "set on copy warning"
    inDf['qstart'], inDf['qend'] = inDf[['qstart','qend']].min(axis=1), inDf[['qstart','qend']].max(axis=1)
    inDf['percmatch']     = (inDf['length'] / inDf['slen']*100)
    inDf['abs percmatch'] = 100 - abs(100 - inDf['percmatch'])#eg changes 102.1->97.9
    inDf['pi_permatch']   = (inDf["pident"] * inDf["abs percmatch"])/100
    inDf['score']         = (inDf['pi_permatch']/100) * inDf["length"]
    inDf['fragment']      = inDf["percmatch"] < 95

    if is_linear == False:
        inDf['qlen']      = (inDf['qlen']/2).astype('int')

    #applies a bonus for anything that is a 100% match to database
    #heurestic! change value maybe
    bonus = 1
    inDf.loc[inDf['pi_permatch']==100, "score"] = inDf.loc[inDf['pi_permatch']==100,'score'] * bonus
    if task == "BLAST": #gives edge to nuc database 
        inDf['score']   = inDf['score'] * 1.1

    wiggleSize = 0.15 #this is the percent "trimmed" on either end eg 0.1 == 90%
    inDf['wiggle'] = (inDf['length'] * wiggleSize).astype(int)
    inDf['wstart'] =  inDf['qstart'] + inDf['wiggle']
    inDf['wend']   =  inDf['qend']   - inDf['wiggle']

    return inDf

def clean(inDf):
    #subtracts a full plasLen if longer than tot length
    inDf['qstart'] = np.where(inDf['qstart'] >= inDf['qlen'], inDf['qstart'] - inDf['qlen'], inDf['qstart'])
    inDf['qend']   = np.where(inDf['qend']   >= inDf['qlen'], inDf['qend']   - inDf['qlen'], inDf['qend'])

    inDf['wstart'] = np.where(inDf['wstart'] >= inDf['qlen'], inDf['wstart'] - inDf['qlen'], inDf['wstart'])
    inDf['wend']   = np.where(inDf['wend']   >= inDf['qlen'], inDf['wend']   - inDf['qlen'], inDf['wend'])

    inDf=inDf.drop_duplicates()
    inDf=inDf.reset_index(drop=True)

    #st.write("raw", inDf)

    #I *think* this has to go before seqspace calcs, but I dont remember the logic
    #inDf=calc_level(inDf)

    #create a conceptual sequence space
    seqSpace=[]
    end    = int(inDf['qlen'][0])

    # for some reason some int columns are behaving as floats -- this converts them
    inDf = inDf.apply(pd.to_numeric, errors='ignore', downcast = "integer")

    for i in inDf.index:
        #end    = inDf['qlen'][0]
        wstart = inDf.loc[i]['wstart'] #changed from qstart
        wend   = inDf.loc[i]['wend']   #changed from qend

        sseqid = [inDf.loc[i]['sseqid']]

        if wend < wstart: # if hit crosses ori
            left   = (wend + 1)          * [1]
            center = (wstart - wend - 1) * [0]
            right  = (end  - wstart + 0) * [1]
        else: # if normal
            left   =  wstart             * [0]
            center = (wend - wstart + 1) * [1]
            right  = (end  - wend   - 1) * [0]

        seqSpace.append(sseqid+left+center+right) #index, not append

    seqSpace=pd.DataFrame(seqSpace,columns=['sseqid'] + list(range(0, end)))
    seqSpace=seqSpace.set_index([seqSpace.index, 'sseqid'])

    #filter through overlaps in sequence space
    toDrop=set()
    for i in range(len(seqSpace)):

        if seqSpace.iloc[i].name in toDrop:
            continue #need to test speed

        end    = inDf['qlen'][0] #redundant, but more readable
        qstart = inDf.loc[seqSpace.iloc[i].name[0]]['qstart']
        qend   = inDf.loc[seqSpace.iloc[i].name[0]]['qend']

        #columnSlice=seqSpace.columns[(seqSpace.iloc[i]==1)] #only columns of hit
        if qstart < qend:
            columnSlice = list(range(qstart, qend + 1))
        else:
            columnSlice = list(range(0,qend + 1)) + list(range(qstart, end))

        rowSlice = seqSpace[columnSlice].any(1) #only the rows that are in the columns of hit
        toDrop   = toDrop | set(seqSpace[rowSlice].loc[i+1:].index) #add the indexs below the current to the drop-set

    ####### For keeping 100% matches
    # keep = inDf[inDf['pi_permatch']==100]
    # keep = set(zip(keep.index, keep['sseqid']))
    # st.write(keep)
    # toDrop = toDrop - keep

    seqSpace = seqSpace.drop(toDrop)
    inDf = inDf.loc[seqSpace.index.get_level_values(0)] #needs shared index labels to work
    inDf = inDf.reset_index(drop=True)

    inDf = calc_level(inDf)

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

def get_gbk(inDf,inSeq, is_linear, record = None):
    #this could be passed a more annotated df
    inDf=inDf.reset_index(drop=True)

    #adds a FeatureLocation object so it can be used in gbk construction
    inDf['feat loc']=inDf.apply(FeatureLocation_smart, axis=1)

    #make a record if one is not provided
    if record is None:
        record = SeqRecord(seq=Seq(inSeq),name='pLannotate')

    if is_linear:
        record.annotations["topology"] = "linear"
    else:
        record.annotations["topology"] = "circular"

    inDf['Type'] = inDf['Type'].str.replace("origin of replication", "rep_origin")
    for index in inDf.index:
        record.features.append(SeqFeature(
            inDf.loc[index]['feat loc'],
            type = inDf.loc[index]["Type"], #maybe change 'Type'
            qualifiers = {
                "note": "pLannotate",
                "label": inDf.loc[index]["Feature"],
                "database":inDf.loc[index]["db"],
                "identity": inDf.loc[index]["pident"],
                "match_length": inDf.loc[index]["percmatch"],
                "fragment": inDf.loc[index]["fragment"],
                "other": inDf.loc[index]["Type"]})) #maybe change 'Type'

    #converts gbk into straight text
    outfileloc=NamedTemporaryFile()
    with open(outfileloc.name, "w") as handle:
        record.annotations["molecule_type"] = "DNA"
        SeqIO.write(record, handle, "genbank")
    with open(outfileloc.name) as handle:
        record=handle.read()
    outfileloc.close()

    return record

def annotate(inSeq, blast_database, linear = False):

    progressBar = st.progress(0)
    progressBar.progress(5)

    #This catches errors in sequence via Biopython
    fileloc = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(inSeq),name="pLannotate",annotations={"molecule_type": "DNA"}), fileloc.name, 'fasta')
    record=list(SeqIO.parse(fileloc.name, "fasta"))
    fileloc.close()

    record=record[0]

    # doubles sequence for origin crossing hits
    if linear == False:
        query = str(record.seq) + str(record.seq)
    elif linear == True:
        query = str(record.seq)
    else:
        progressBar.empty()
        st.error("error")
        return pd.DataFrame()

    #addgene BLAST
    database = os.path.join(blast_database, "addgene_collected_features_test_20-12-11")
    nucs = BLAST(seq=query, wordsize=12, db=database, task = "BLAST")
    nucs = calculate(nucs, task = "BLAST", is_linear = linear)
    nucs['db'] = "addgene"

    progressBar.progress(25)

    #orfs = find_orfs(query, linear)
    database=" ".join(os.path.join(blast_database, x) for x in ("Rfam.clanin", "Rfam.cm"))
    rnas = BLAST(seq=query, wordsize=12, db=database, task = "infernal")
    rnas['qlen'] = len(query)
    rnas = calculate(rnas, task = "infernal", is_linear = linear)
    rnas['db'] = "infernal"

    progressBar.progress(55)

    #swissprot DIAMOND search
    database=os.path.join(blast_database, "trimmed_swissprot.dmnd")
    prots = BLAST(seq=query,wordsize=12, db=database, task="DIAMOND")
    prots = calculate(prots, task = "DIAMOND", is_linear = linear) #calc not explicit
    prots['db'] = "swissprot"

    progressBar.progress(75)

    #fpbase DIAMOND search
    database=os.path.join(blast_database, "fpbase.dmnd")
    fluors = BLAST(seq=query,wordsize=12, db=database, task="DIAMOND")
    fluors = calculate(fluors, task = "DIAMOND", is_linear =  linear) #calc not explicit
    fluors['db'] = "fpbase"

    progressBar.progress(90)

    #aggregates all dfs together and sorts
    blastDf = nucs.append(prots)
    blastDf = blastDf.append(fluors)
    blastDf = blastDf.append(rnas)
    #blastDf = blastDf.append(orfs)
    blastDf = blastDf.sort_values(by=["score","length","percmatch"], ascending=[False, False, False])

    if blastDf.empty: #if no hits are found
        progressBar.empty()
        return blastDf

    blastDf = clean(blastDf)

    if blastDf.empty: #if no hits are found
        progressBar.empty()
        return blastDf

    progressBar.empty()

    blastDf['blastDf'] = blastDf['qend'] + 1 #corrects position for gbk

    #blastDf = blastDf.append(orfs)

    return blastDf

