import numpy as np
import pandas as pd

def parse_infernal(file_loc):

    with open(file_loc) as file_handle:
        lines = file_handle.readlines()

    #find position of columns using "---" field
    #create list of len == 2 tuples
    col_widths = [len(ele)+1 for ele in lines[1].split()]
    ends = list(np.cumsum(col_widths))
    #changes longest data line in file
    #ends[-1] = max(map(len,[line for line in lines if line[0] != "#"]))
    ends[-1] += 100 #not super elegant -- just adds 100 to capture full output
    starts = [0] + ends[:-1]
    col_pos = list(zip(starts,ends))

    #extract column names using above positions
    col_names = []
    for ele in col_pos:
        col_names.append(lines[0][ele[0]:ele[1]].strip())

    try:
        infernal = pd.read_fwf(file_loc, comment="#", colspecs = col_pos, header = None)
        infernal.columns = col_names
    except pd.errors.EmptyDataError:
        infernal = pd.DataFrame(columns = col_names)

    columns = ['#idx','target name', 'accession','clan name','seq from', 'seq to','mdl from','mdl to','strand', 'score', 'E-value','description of target']
    infernal = infernal[columns]
    infernal = infernal.loc[:,~infernal.columns.duplicated()]
    replacements = {"#idx":"sseqid","seq from":"qstart","seq to":"qend",'mdl from':"sstart",'mdl to':"send","E-value":"evalue","strand":"sframe"}
    infernal = infernal.rename(columns=replacements)
    infernal["accession"] = infernal["accession"].str.replace("-"," ")
    infernal["clan name"] = infernal["clan name"].str.replace("-"," ")

    infernal = infernal.rename(columns = {"target name":"Feature","description of target":"Description"})
    infernal['Feature'] = infernal['Feature'].str.replace("_"," ")
    infernal['Description'] = "Accession: " + infernal['accession'] + " - " + infernal['Description']
    
    infernal = infernal.apply(pd.to_numeric, errors='ignore', downcast = "integer")
    
    infernal['qseq'] = ""
    infernal["sframe"] = infernal["sframe"].replace(["-","+"], [-1,1])
    infernal['qstart'] = infernal['qstart']-1
    infernal['qend']   = infernal['qend']-1
    infernal['length'] = abs(infernal['qend']-infernal['qstart'])+1
    infernal['slen'] = abs(infernal['send']-infernal['sstart'])+1
    infernal['pident'] = 100

    #clan name currently not used
    infernal = infernal.drop(columns=['accession','clan name'])
    
    return infernal
