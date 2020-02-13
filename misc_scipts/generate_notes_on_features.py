from Bio import SeqIO
import glob
import pandas as pd

c=[]
for infile_loc in glob.glob('/Users/mattmcguffie/database/all_snapgene_features/*.gb'):
    with open(infile_loc,'r') as file_handle:
        record_dict = SeqIO.to_dict(SeqIO.parse(file_handle, 'gb'))

    fileName=infile_loc.split("/")[-1].split(".gb")[0]
    gbkFile = record_dict[list(record_dict.keys())[0]]

    featureType=gbkFile.features[1].type.replace("_"," ")
    if featureType == "rep origin":
        featureType = "origin of replication"

    try: l = gbkFile.features[1].qualifiers['label'][0]
    except: l = None
    try: n = gbkFile.features[1].qualifiers['note'][0]
    except: n = None
    try: p = gbkFile.features[1].qualifiers['product'][0]
    except: p = None

    c.append((fileName,l,featureType,n,p))

df=pd.DataFrame(c,columns=["file_name","label","type","note","product"])

desc_df=[]

for i in range(len(df)):
    fileName=df.iloc[i]['file_name']
    featureType=df.iloc[i]['type']
    note=df.iloc[i]['note']
    prod=df.iloc[i]['product']
    label=df.iloc[i]['label']
    if note and prod:
        desc=f"{prod.capitalize()} — {note}"
    elif note:
        desc=f"{note.capitalize()}"
    elif prod:
        desc=f"{prod.capitalize()}"
    else:
        desc=featureType

    if desc[-1]!=".":
        desc+="."

    desc_df.append((fileName,label,featureType,desc))

desc_df=pd.DataFrame(desc_df,columns=["file","Feature","Type","Description"]).set_index("file",drop=True)

desc_df.to_csv("./feature_notes.csv",sep="\t")
