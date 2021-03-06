import streamlit as st
import numpy as np
import base64
from Bio import SeqIO
from tempfile import NamedTemporaryFile
import pandas as pd
from annotate import annotate, get_gbk
from bokeh_plot import get_bokeh
import glob
from BLAST_hit_details import details
import io
import sys

#sys.tracebacklimit = 0 #removes traceback so code is not shown during errors

hide_streamlit_style = """
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style>
"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True) 

st.image("./images/pLannotate.png",use_column_width=False, width=500)

#markdown hack to remove full screen icon from pLannotate logo
hide_full_screen = '''
<style>
.element-container:nth-child(2) button{visibility: hidden;}
</style>
'''
st.markdown(hide_full_screen, unsafe_allow_html=True) 

st.subheader('Annotate your engineered plasmids')
sidebar = st.sidebar.empty()
sidebar.markdown('''
        **<a href="https://en.wikipedia.org/wiki/Plasmid" target="_blank">Plasmids</a>** are ubiquitous in many fields of biology.

        Engineered plasmids generally have long and circuitous cloning histories, meaning annotations are forgotten, and often contain cryptic genes and gene fragments leftover from cloning.

        **<font color="#f9a557">pLannotate</font>** re-annotates engineered plasmids and shows you where the fragments are.''',unsafe_allow_html=True)

inSeq=""
maxPlasSize = 50000
IUPAC= 'GATCRYWSMKHBVDNgatcrywsmkhbvdn'
option = st.radio(
    'Choose method of submitting sequence:',
    ["Upload a file (.fa or .fasta)", "Enter a sequence","Example"])

if option == "Upload a file (.fa or .fasta)":

    #markdown css hack to remove fullscreen -- fickle because it is hardcoded
    nth_child_num = 13

    uploaded_file = st.file_uploader("Choose a file:", type=['fa',"fasta"])

    if uploaded_file is not None:
        text_io = io.TextIOWrapper(uploaded_file,encoding='UTF-8')

        st.success("File uploaded.")
        
        #This catches errors on file uploads via Biopython
        fileloc = NamedTemporaryFile()
        record = list(SeqIO.parse(text_io, "fasta"))
        record[0].annotations["molecule_type"] = "DNA"
        SeqIO.write(record, fileloc.name, 'fasta')
        record = list(SeqIO.parse(fileloc.name, "fasta"))
        fileloc.close()

        if len(record)!=1:
            error = 'FASTA file contains many entries --> please submit a single FASTA file.'
            raise ValueError(error)
        
        inSeq = str(record[0].seq)

elif option == "Enter a sequence":

    #markdown css hack to remove fullscreen -- fickle because it is hardcoded
    nth_child_num = 12

    inSeq = st.text_area('Input sequence here:',max_chars = maxPlasSize)
    inSeq = inSeq.replace("\n","")

elif option == "Example":
    
    #markdown css hack to remove fullscreen -- fickle because it is hardcoded
    nth_child_num = 12
    
    fastas=[]
    for infile_loc in glob.glob('./fastas/*.fa'):
        fastas.append(infile_loc.split("/")[-1].split(".fa")[0])
    exampleFile = st.radio("Choose example file:",fastas)
    inSeq = str(list(SeqIO.parse(f"./fastas/{exampleFile}.fa", "fasta"))[0].seq)

if inSeq:

    if not set(inSeq).issubset(IUPAC):
        error = f'Sequence contains invalid characters -- must be ATCG and/or valid IUPAC nucleotide ambiguity code'
        raise ValueError(error)

    if len(inSeq) > maxPlasSize:
        error = f'Are you sure this is an engineered plasmid? Entry size is too large -- must be {maxPlasSize} bases or less.'
        raise ValueError(error)

    with st.spinner("Annotating..."):
        recordDf = annotate(inSeq)

        if recordDf.empty:
            st.error("No annotations found.")
        else:
            with open("./FAQ.md") as fh:
                faq = fh.read()
            sidebar.markdown(faq)
            recordDf = details(recordDf)

            st.markdown("---")
            st.header('Results:')
            
            st.write("Hover mouse for info, click and drag to pan, scroll wheel to zoom")
            st.bokeh_chart(get_bokeh(recordDf), use_container_width=False)
            
            #markdown hack to remove full screen icon from bokeh plot (which is distorted)
            hide_full_screen = f'''
            <style>
            .element-container:nth-child({nth_child_num}) button{{visibility: hidden;}}
            </style>
            '''
            st.markdown(hide_full_screen, unsafe_allow_html=True) 

            st.header("Download Annotations:")

            #creates a procedurally-gen name for file based on seq 
            filename = str(abs(hash(inSeq)))[:6]

            #write and encode gbk for dl
            gbk=get_gbk(recordDf,inSeq)
            b64 = base64.b64encode(gbk.encode()).decode()
            gbk_dl = f'<a href="data:text/plain;base64,{b64}" download="{filename}.gbk"> download {filename}.gbk</a>'
            st.markdown(gbk_dl, unsafe_allow_html=True)

            #encode csv for dl
            columns = ['qstart', 'qend', 'sframe', 'pident', 'slen', 'sseq', 'length', 'uniprot', 'abs percmatch', 'fragment', 'db', 'Feature', 'Type', 'Description']
            replacements = {'qstart':'start location', 'qend':'end location', 'sframe':'strand', 'pident':'percent identity', 'slen':'full length of feature in db', 'sseq':'full sequence of feature in db', 'length':'length of found feature', 'uniprot':'uniprot ID', 'abs percmatch':'percent match length', 'db':'database'}
            cleaned = recordDf[columns]
            cleaned = cleaned.rename(columns=replacements)
            csv = cleaned.to_csv(index=False)
            b64 = base64.b64encode(csv.encode()).decode()
            csv_dl = f'<a href="data:text/plain;base64,{b64}" download="{filename}.csv"> download {filename}.csv</a>'
            st.markdown(csv_dl, unsafe_allow_html=True)

            st.markdown("---")

            #prints table of features
            st.header("Features")
            displayColumns = ['Feature','percent identity','percent match length','Description']
            markdown = cleaned[displayColumns].copy()
            numericCols = ['percent identity', 'percent match length']
            markdown[numericCols] = np.round(markdown[numericCols], 1)
            markdown[numericCols] = markdown[numericCols].astype(str) + "%"
            markdown = markdown.set_index("Feature",drop=True)
            st.markdown(markdown.drop_duplicates().to_markdown())