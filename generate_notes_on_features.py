from Bio import SeqIO
import glob
import pandas as pd

c=[]
for infile_loc in glob.glob('/Users/mattmcguffie/database/all_snapgene_features/*.gb'):
    with open(infile_loc,'r') as file_handle:
        record_dict = SeqIO.to_dict(SeqIO.parse(file_handle, 'gb'))
        
    fileName=infile_loc.split("/")[-1]
    
    gbkFile = record_dict[list(record_dict.keys())[0]]
    try: l = gbkFile.features[1].qualifiers['label'][0]
    except: l = None
    try: n = gbkFile.features[1].qualifiers['note'][0]
    except: n = None
    try: p = gbkFile.features[1].qualifiers['product'][0]
    except: p = None
    
    c.append((fileName,l,n,p))
    
df=pd.DataFrame(c,columns=["file_name","label","note","product"]).set_index(["file_name"],drop=True)
df