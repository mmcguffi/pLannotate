import subprocess
from tempfile import NamedTemporaryFile

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd
import streamlit as st

import plannotate.resources as rsc
from plannotate.infernal import parse_infernal

log = NamedTemporaryFile()

def BLAST(seq, db):
    task = db['method']
    parameters = db['parameters']
    db_loc = db['db_loc']
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    if task == "blastn":
        flags = 'qstart qend sseqid sframe pident slen qseq length sstart send qlen evalue'
        subprocess.call( 
            (f'blastn -task blastn-short -query {query.name} -out {tmp.name} ' 
             f'-db {db_loc} {parameters} -outfmt "6 {flags}" >> {log.name} 2>&1'),
            shell=True)

    elif task == "diamond":
        flags = 'qstart qend sseqid pident slen qseq length sstart send qlen evalue'
        subprocess.call(f'diamond blastx -d {db_loc} -q {query.name} -o {tmp.name} '
                        f'{parameters} --outfmt 6 {flags} >> {log.name} 2>&1',shell=True)

    elif task == "infernal":
        flags = "--cut_ga --rfam --noali --nohmmonly --fmt 2" 
        cmd = f"cmscan {flags} --tblout {tmp.name} --clanin {db_loc} {query.name} >> {log.name} 2>&1"
        subprocess.call(cmd, shell=True)
        inDf = parse_infernal(tmp.name)
        
        inDf['qlen'] = len(seq)
        
        #manually gets DNA sequence from seq(x2)
        if not inDf.empty:
            inDf['qseq'] = inDf.apply(lambda x: (seq)[x['qstart']:x['qend']+1].upper(), axis=1)
        
        tmp.close()
        query.close()

        return inDf

    with open(tmp.name, "r") as file_handle:  #opens BLAST file
        align = file_handle.readlines()

    tmp.close()
    query.close()

    inDf = pd.DataFrame([ele.split() for ele in align],columns=flags.split())
    inDf = inDf.apply(pd.to_numeric, errors='ignore')

    if task == "diamond":
        try:
            inDf['sseqid'] = inDf['sseqid'].str.split("|", n=2, expand=True)[1]
        except (ValueError, KeyError):
            pass
        inDf['sframe'] = (inDf['qstart']<inDf['qend']).astype(int).replace(0,-1)
        inDf['slen']   = inDf['slen'] * 3
        inDf['length'] = abs(inDf['qend']-inDf['qstart'])+1

    return inDf

def calculate(inDf, is_linear):

    inDf['qstart'] = inDf['qstart']-1
    inDf['qend']   = inDf['qend']-1

    inDf['qstart'], inDf['qend'] = inDf[['qstart','qend']].min(axis=1), inDf[['qstart','qend']].max(axis=1)
    inDf['percmatch']     = (inDf['length'] / inDf['slen']*100)
    inDf['abs percmatch'] = 100 - abs(100 - inDf['percmatch'])#eg changes 102.1->97.9
    inDf['pi_permatch']   = (inDf["pident"] * inDf["abs percmatch"])/100
    inDf['score']         = (inDf['pi_permatch']/100) * inDf["length"]
    
    # score adjustment heuristic
    # higher priority == less score deduction
    # each prirority num increase decreases score by 1/2
    # eg: priority 1 == 1 | priority 2 == 1/2 | priority 3 == 1/4 | etc
    inDf['score']  = inDf['score'] * (2**(-1 * inDf['priority'].astype(float)) * 2)

    if is_linear == False:
        inDf['qlen']      = (inDf['qlen']/2).astype('int')

    # applies a bonus for anything that is a 100% match to database
    # heurestic! bonus depends on priority
    bonus = (1/inDf['priority']) * 10
    inDf.loc[inDf['pi_permatch']==100, "score"] = inDf.loc[inDf['pi_permatch']==100,'score'] * bonus

    wiggleSize = 0.15 #this is the percent "trimmed" on either end eg 0.1 == 90%
    inDf['wiggle'] = (inDf['length'] * wiggleSize).astype(int)
    inDf['wstart'] =  inDf['qstart'] + inDf['wiggle']
    inDf['wend']   =  inDf['qend']   - inDf['wiggle']

    return inDf

def clean(inDf):
    #subtracts a full plasLen if longer than tot length
    inDf['qstart_dup'] = inDf['qstart']
    inDf['qend_dup']   = inDf['qend']
    inDf['qstart'] = np.where(inDf['qstart'] >= inDf['qlen'], inDf['qstart'] - inDf['qlen'], inDf['qstart'])
    inDf['qend']   = np.where(inDf['qend']   >= inDf['qlen'], inDf['qend']   - inDf['qlen'], inDf['qend'])

    inDf['wstart'] = np.where(inDf['wstart'] >= inDf['qlen'], inDf['wstart'] - inDf['qlen'], inDf['wstart'])
    inDf['wend']   = np.where(inDf['wend']   >= inDf['qlen'], inDf['wend']   - inDf['qlen'], inDf['wend'])

    # these are manually-curated (garbage) hits that overlap with common features
    problem_hits = ['P03851','P03845','ISS','P03846']
    inDf = inDf.loc[~inDf['sseqid'].isin(problem_hits)]
    
    # filter for evalue less than 1 (should only affect SnapGene db?)
    inDf = inDf.loc[inDf['evalue'] < 1]
    
    # drop poor matches that are very small fragments
    # usually an artifact from wonky SnapGene features that are composite features
    inDf = inDf.loc[inDf['pi_permatch'] > 3]

    inDf=inDf.drop_duplicates()
    inDf=inDf.reset_index(drop=True)
    
    if inDf.empty:
        return inDf
        
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
            left   = (wend + 1)          * [inDf.loc[i]['kind']]
            center = (wstart - wend - 1) * [None]
            right  = (end  - wstart + 0) * [inDf.loc[i]['kind']]
        else: # if normal
            left   =  wstart             * [None]
            center = (wend - wstart + 1) * [inDf.loc[i]['kind']]
            right  = (end  - wend   - 1) * [None]

        seqSpace.append(sseqid+left+center+right) #index, not append

    seqSpace=pd.DataFrame(seqSpace,columns=['sseqid'] + list(range(0, end)))
    seqSpace=seqSpace.set_index([seqSpace.index, 'sseqid']) #multi-indexed
    #filter through overlaps in sequence space
    toDrop=set()
    for i in range(len(seqSpace)):

        if seqSpace.iloc[i].name in toDrop:
            continue #need to test speed

        end    = inDf['qlen'][0] #redundant, but more readable
        qstart = inDf.loc[seqSpace.iloc[i].name[0]]['qstart']
        qend   = inDf.loc[seqSpace.iloc[i].name[0]]['qend']
        kind   = inDf.loc[seqSpace.iloc[i].name[0]]['kind']

        #columnSlice=seqSpace.columns[(seqSpace.iloc[i]==1)] #only columns of hit
        if qstart < qend:
            columnSlice = list(range(qstart+1, qend + 1))
        else:
            columnSlice = list(range(0,qend + 1)) + list(range(qstart, end))
        
        rowSlice = (seqSpace[columnSlice] == kind).any(1) #only the rows that are in the columns of hit
        toDrop   = toDrop | set(seqSpace[rowSlice].loc[i+1:].index) #add the indexs below the current to the drop-set

    seqSpace = seqSpace.drop(toDrop)
    inDf = inDf.loc[seqSpace.index.get_level_values(0)] #needs shared index labels to work
    inDf = inDf.reset_index(drop=True)
    # may need to run this with df that "passes" the origin

    return inDf


def get_details(inDf, yaml_file_loc):
    
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
    databases = rsc.get_yaml(yaml_file_loc)
    
    assert len(set(inDf['db'].to_list())) == 1, "All hits must be from the same database"
    database_name = inDf['db'].to_list()[0]
    
    
    database = databases[database_name]
    
    sseqids = inDf.loc[inDf['db'] == database_name]['sseqid'].tolist()
    sseqids = [_ for _ in sseqids if _] #removes blank edgecases
    
    # this manually exctracts "3xHA" from "pdb|3xHA|"
    # probably other instances of this issue, cannot track down source of this issue
    # pretty hacky, but it works
    problem_name = r"pdb\|(.*)\|"
    inDf['sseqid'] = inDf['sseqid'].str.replace(problem_name, r"\1", regex=True)

    db_details = database['details']

    if db_details['location'] == 'None':
        #if no file is passed, data should already be in dataframe
        feat_desc = inDf.loc[inDf['db'] == database_name][['sseqid','Feature','Description']]
        
    else:
        if db_details['location'] == "Default":
            details_file_loc = rsc.get_details(database_name) + ".csv"
        else: #if a file path is passed, use that
            details_file_loc = db_details['location']

        #if the description file is compressed
        if db_details['compressed'] == True:
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
            # if 'existence level' is not found, 
            # 0 is returned as the location
            # meaning 15 and 16 are the default values
            # this sets a baseline priority of `1` if nothing is found
                if s == 15 and e == 16:
                    return 0
                else:
                    return int(d[s:e]) - 1

            # extract the level from the description
            feat_desc['priority_mod'] = [calc_priority_mod(d,s,e) for d, s, e in zip(feat_desc["Description"], feat_desc["s"], feat_desc["e"])]
            feat_desc = feat_desc.drop(columns=['s','e'])

    #try to see if a default type was passed
    if db_details['default_type'] != 'None':
        feat_desc['Type'] = db_details['default_type']
    else:
        pass 
    
    return feat_desc


@st.cache(hash_funcs={pd.DataFrame: lambda _: None}, suppress_st_warning=True, max_entries = 10, show_spinner=False)
def get_raw_hits(query, linear, yaml_file_loc):
     
    progressBar = st.progress(0)
    progress_amt = 5
    progressBar.progress(progress_amt)

    databases = rsc.get_yaml(yaml_file_loc)
    increment = int(90 / len(databases))
    
    raw_hits = []
    for database_name in databases:
        database = databases[database_name]
        hits = BLAST(seq = query, db = database)
        
        hits['db'] = database_name
        hits['sseqid'] = hits['sseqid'].astype(str)
        
        if hits.empty:
            continue
        
        feat_descriptions = get_details(hits, yaml_file_loc)
        # `suffixes = ('_x', None)` means the descriptions for Rfam will be copied,
        # the original descriptions will be appeneded with `_x` and can be ignored
        # the Rfam descriptions are in the original df due to the quirks of how the details
        # are stored, so this is a work around. Possibly condsider dropping the `_x`` column
        hits = hits.merge(feat_descriptions, on='sseqid', how='left', suffixes = ('_x', None))
        hits = hits[hits.columns.drop(list(hits.filter(regex='_x')))]
        
        #removes primer binding site annotations
        hits = hits.loc[hits['Type'] != 'primer_bind']
        
        hits['priority'] = database['priority']
        try:
            hits['priority'] = hits['priority'] + hits['priority_mod']
            hits = hits.drop('priority_mod', axis=1)
        except KeyError:
            pass
        hits = calculate(hits, is_linear = linear)
                
        raw_hits.append(hits)
        
        progress_amt += increment
        progressBar.progress(progress_amt)
        
    if len(raw_hits) == 0:
        return pd.DataFrame()
    
    blastDf = pd.concat(raw_hits)
    
    blastDf = blastDf.sort_values(by=["score","length","percmatch"], ascending=[False, False, False])

    progressBar.empty()
    
    return blastDf

def annotate(inSeq, yaml_file = rsc.get_yaml_path(), linear = False, is_detailed = False):

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
        st.error("error")
        return pd.DataFrame()
    
    blastDf = get_raw_hits(query, linear, yaml_file)
    
    if blastDf.empty: #if no hits are found
        return blastDf
    
    #this has to re-parse the yaml, so not an elegant solution
    if is_detailed == True:
        blastDf['kind'] = blastDf['Type']
    else:
        blastDf['kind'] = 1
    
    blastDf = clean(blastDf)
    
    if blastDf.empty: #if no hits are found
        return blastDf
    
    def is_fragment(feature):
        if feature['Type'] == "CDS":
            if feature['pi_permatch'] == 100:
                return False
            elif ((feature['length'] % 3) == 0) & (feature["percmatch"] > 95):
                return False
            else:
                return True
        elif feature['Type'] != "CDS":
            if feature['percmatch'] < 95:
                return True
            else:
                return False
        else:
            st.error("Fragment error.")
    blastDf['fragment'] =  blastDf.apply(is_fragment, axis=1)
    
    if blastDf.empty: #if no hits are found
        return blastDf

    blastDf['qend'] = blastDf['qend'] + 1 #corrects position for gbk

    #manually gets DNA sequence from inSeq
    #blastDf['qseq'] = inSeq #adds the sequence to the df
    #blastDf['qseq'] = blastDf.apply(lambda x: x['qseq'][x['qstart']:x['qend']+1], axis=1)
    blastDf['qseq'] = blastDf.apply(lambda x: str(Seq(x['qseq']).reverse_complement()) if x['sframe'] == -1 else x['qseq'], axis=1)

    global log
    log.close()

    return blastDf