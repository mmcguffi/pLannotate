from Bio.Seq import Seq
import orfipy_core 
import pandas as pd
import streamlit as st

def find_orfs(plasmid_seq, is_linear):
    size = 300
    start_codons = ["ATG"]
    orfs = []
    try:
        hits = orfipy_core.orfs(str(plasmid_seq).upper(),minlen=size,starts=start_codons,include_stop=True)
    except IndexError:
        return pd.DataFrame()
    for orf in hits:
        qlen = len(plasmid_seq)
        if is_linear == False:
            qlen = qlen/2
        sseqid = orf[3].split(";")[0]
        start = orf[0]
        end = orf[1]
        strand = orf[2]
        if strand == '-':
            sframe = -1
        else:
            sframe = 1
        orfs.append((qlen,sseqid,start,end,sframe,"ORF"))
        # seq = plas[start:end]
        # if strand == '-':
        #     seq = seq.reverse_complement()
    orfs = pd.DataFrame(orfs,columns=["qlen","sseqid","qstart","qend","sframe","db"])
    orfs[["wstart","wend"]] = orfs[['qstart','qend']]
    #st.write(orfs)
    orfs = orfs.apply(pd.to_numeric, errors='ignore', downcast = "integer")
    return orfs
