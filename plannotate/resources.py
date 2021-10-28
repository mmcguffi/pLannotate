
import os
from tempfile import NamedTemporaryFile

import pkg_resources
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

valid_genbank_exts = ['.gbk', '.gb', '.gbf', '.gbff']
valid_fasta_exts = ['.fa', '.fasta']
maxPlasSize = 50000


def get_resource(group, name):
    return pkg_resources.resource_filename(__package__, f"data/{group}/{name}")


def get_image(name):
    return get_resource("images", name)


def get_template(name):
    return get_resource("templates", name)


def get_example_fastas():
    return get_resource("fastas", "")


def get_yaml(blast_database_loc):
    return parse_yaml(get_resource("data", "databases.yml"),blast_database_loc)

def get_details(name):
    return get_resource("data", name)


def get_name_ext(file_loc):
    base = os.path.basename(file_loc)
    name = os.path.splitext(base)[0]
    ext = os.path.splitext(base)[1]
    return name,ext


def validate_file(file, ext):
    if ext in valid_fasta_exts:
        #This catches errors on file uploads via Biopython
        temp_fileloc = NamedTemporaryFile()
        record = list(SeqIO.parse(file, "fasta"))
        try:
            record[0].annotations["molecule_type"] = "DNA"
        except IndexError:
            error = "Malformed fasta file --> please submit a fasta file in standard format"
            raise ValueError(error)
        SeqIO.write(record, temp_fileloc.name, 'fasta')
        record = list(SeqIO.parse(temp_fileloc.name, "fasta"))
        temp_fileloc.close()

        if len(record)!=1:
            error = 'FASTA file contains many entries --> please submit a single FASTA file.'
            raise ValueError(error)

    elif ext in valid_genbank_exts:
        temp_fileloc = NamedTemporaryFile()
        try:
            record = list(SeqIO.parse(file, "gb"))[0]
        except IndexError:
            error = "Malformed Genbank file --> please submit a Genbank file in standard format"
            raise ValueError(error)
        # submitted_gbk = record # for combining -- not current imlementated
        SeqIO.write(record, temp_fileloc.name, 'fasta')
        record = list(SeqIO.parse(temp_fileloc.name, "fasta"))
        temp_fileloc.close()
    
    else:
        error = 'must be a FASTA or GenBank file'
        raise ValueError(error)

    if len(record)!=1:
        error = 'FASTA file contains many entries --> please submit a single FASTA file.'
        raise ValueError(error)
        
    inSeq = str(record[0].seq)

    validate_sequence(inSeq)

    return inSeq


def validate_sequence(inSeq):
    IUPAC= 'GATCRYWSMKHBVDNgatcrywsmkhbvdn'
    if not set(inSeq).issubset(IUPAC):
        error = f'Sequence contains invalid characters -- must be ATCG and/or valid IUPAC nucleotide ambiguity code'
        raise ValueError(error)

    if len(inSeq) > maxPlasSize:
        error = f'Are you sure this is an engineered plasmid? Entry size is too large -- must be {maxPlasSize} bases or less.'
        raise ValueError(error)


def get_gbk(inDf,inSeq, is_linear, record = None):
    #this could be passed a more annotated df
    inDf=inDf.reset_index(drop=True)

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


def get_clean_csv_df(recordDf):
    # change sseqid to something more legible
    columns = ['sseqid', 'qstart', 'qend', 'sframe', 'pident', 'slen', 'length', 'abs percmatch', 'fragment', 'db', 'Feature', 'Type', 'Description', 'qseq']
    cleaned = recordDf[columns]
    replacements = {'qstart':'start location', 'qend':'end location', 'sframe':'strand', 'pident':'percent identity', 'slen':'full length of feature in db', 'qseq':'sequence', 'length':'length of found feature', 'abs percmatch':'percent match length', 'db':'database'}
    cleaned = cleaned.rename(columns=replacements)
    return cleaned

import yaml

#parse yaml file
# def parse_yaml(file_name):
#     with open(file_name, 'r') as f:
#         dbs = yaml.load(f, Loader = yaml.SafeLoader)

#     for db in dbs.keys():
#         method = dbs[db]['method']
#         try:
#             parameters = " ".join(dbs[db]['parameters'])
#         except KeyError:
#             parameters = ""
#         details = dbs[db]['details']
#         #print(f'{method} {parameters} {details}')
#         return method, parameters, details

#         print()

def parse_yaml(file_name, blast_database_loc):
    with open(file_name, 'r') as f:
        dbs = yaml.load(f, Loader = yaml.SafeLoader)

    #collapes list
    db_list = []
    for db in dbs.keys():
        try:
            parameters = " ".join(dbs[db]['parameters'])
        except KeyError:
            parameters = ""
        dbs[db]['parameters'] = parameters

        if dbs[db]['method'] == 'infernal':
            db_loc=" ".join(os.path.join(blast_database_loc, x) for x in (f"{db}.clanin", f"{db}.cm"))
        else:
            db_loc = os.path.join(blast_database_loc, db)
        #dbs[db]['name'] = db
        dbs[db]['db_loc'] = db_loc
        #db_list.append(dbs[db])
    return dbs