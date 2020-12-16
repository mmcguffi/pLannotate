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
from io import StringIO   

st.image("./images/pLannotate.png",use_column_width=True, width=500)
st.subheader('v0.3')

st.sidebar.markdown('''
        **<a href="https://en.wikipedia.org/wiki/Plasmid" target="_blank">Plasmids</a>** are ubiquitous in many fields of biology.

        Engineered plasmids generally have long and circuitous cloning histories, meaning annotations are forgotten, and often contain hidden junk leftover from cloning.

        **<font color="#f9a557">pLannotate</font>** re-annotates engineered plasmids and shows you where the junk is.''',unsafe_allow_html=True)

inSeq=""

option = st.radio(
    'Choose method of submitting sequence:',
    ["Upload a file (.fa or .fasta)", "Enter a sequence","Example"])

# def get_uploaded_file():
#     #taken from: https://github.com/streamlit/streamlit/issues/2266
#     file_buffer = st.file_uploader("Choose a file:", type=['fa',"fasta"])
#     if file_buffer is None:
#         return None
#     file_buffer.seek(0)
#     uploaded_file = io.TextIOWrapper(file_buffer,encoding='UTF-8')
#     return(uploaded_file)


if option == "Upload a file (.fa or .fasta)":
    #this is to supress a decrecation warning -- not best practice
    #st.set_option('deprecation.showfileUploaderEncoding', False)
    #uploaded_file = st.file_uploader("Choose a file:", type=['fa',"fasta"],encoding='UTF-8')

    # file = get_uploaded_file()
    # if file is not None:
    #     st.success("File uploaded.")
    #     file=file.readlines()
    #     inSeq="".join(file[1:]).strip().replace("\n","").replace("\r","")

    uploaded_file = st.file_uploader("Choose a file:", type=['fa',"fasta"])

    if uploaded_file is not None:
        text_io = io.TextIOWrapper(uploaded_file,encoding='UTF-8')
        #uploaded_file = StringIO(uploaded_file.decode("utf-8"))
        st.success("File uploaded.")
        # file=text_io.readlines()
        # inSeq="".join(file[1:]).strip().replace("\n","").replace("\r","")
        
        #I could just create a seq object/readlines but this catches errors
        fileloc = NamedTemporaryFile()
        record=list(SeqIO.parse(text_io, "fasta"))
        SeqIO.write(record, fileloc.name, 'fasta')
        record=list(SeqIO.parse(fileloc.name, "fasta"))
        fileloc.close()
        inSeq = str(record[0].seq)

elif option == "Enter a sequence":
    inSeq = st.text_area('Input sequence here:')
    inSeq = inSeq.replace("\n","")
elif option == "Example":
    fastas=[]
    for infile_loc in glob.glob('./fastas/*.fa'):
        fastas.append(infile_loc.split("/")[-1].split(".fa")[0])
    exampleFile=st.radio("Choose example file:",fastas)
    inSeq=str(list(SeqIO.parse(f"./fastas/{exampleFile}.fa", "fasta"))[0].seq)
    #st.text_area('Input sequence here:',inSeq)

if inSeq:
    with st.spinner("Annotating..."):
        recordDf = annotate(inSeq)

        # my_bar = st.progress(0)
        # import time
        # for percent_complete in range(100):
        #     time.sleep(0.1)
        #     my_bar.progress(percent_complete + 1)

        if recordDf.empty:
            st.error("No annotations found.")
        else:

            ######
            recordDf = details(recordDf)
            ######

            st.markdown("---")
            st.header('Results:')

            st.bokeh_chart(get_bokeh(recordDf),use_container_width=False)

            st.header("Download Annotations:")

            # filename = st.text_input('Enter custom file name for download:',"pLannotate")
            # if not filename:
            #     filename="pLannotate"
            filename = str(abs(hash(inSeq)))[:6]

            #write and encode gbk for dl
            gbk=get_gbk(recordDf,inSeq)
            b64 = base64.b64encode(gbk.encode()).decode()

            # dlPicLoc='dl_arrow.png'
            # with open("./dl_arrow.png", "rb") as image_file:
            #     encoded_string = base64.b64encode(image_file.read())

            gbk_dl = f'<a href="data:text/plain;base64,{b64}" download="{filename}.gbk"> download {filename}.gbk</a>'
            st.markdown(gbk_dl, unsafe_allow_html=True)

            #encode csv for dl
            csv = recordDf.to_csv(index=False)
            b64 = base64.b64encode(csv.encode()).decode()
            #csv_dl = f'<a href="data:text/plain;base64,{b64}" download="{filename}.csv"> ![dl]({dlPicLoc} "download .csv") download {filename}.csv</a>'
            csv_dl = f'<a href="data:text/plain;base64,{b64}" download="{filename}.csv"> download {filename}.csv</a>'
            st.markdown(csv_dl, unsafe_allow_html=True)

            st.markdown("---")

            #prints table of features
            frag=recordDf[recordDf['fragment']==True]
            full=recordDf[recordDf['fragment']==False]
            st.header("Features")
            st.markdown(full[['Feature','db','Description']].set_index("Feature",drop=True).to_markdown())
            st.markdown("---")
            if not frag.empty:
                st.header("Possibly Fragmented Features")
                st.markdown(frag[['Feature','db','Description']].set_index("Feature",drop=True).to_markdown())