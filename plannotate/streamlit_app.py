import base64
import glob
import io
import os
import sys

import numpy as np
import streamlit as st
from Bio import SeqIO

import plannotate.resources as rsc
from plannotate.annotate import annotate
from plannotate.BLAST_hit_details import details
from plannotate.bokeh_plot import get_bokeh
from plannotate import __version__ as plannotate_version


def run_streamlit(args):
    
    sidebar, cite_fund, images = setup_page()

    inSeq=""

    # markdown css hack to remove fullscreen
    # fickle because it is hardcoded and can
    # change with streamlit versions updates
    nth_child_num = 14
    
    upload_option  = "Upload a file (FASTA or GenBank)"
    enter_option   = "Enter a sequence"
    example_option = "Example"
    
    option = st.radio(
        'Choose method of submitting sequence:',
        [upload_option, enter_option, example_option]
        )
    
    if option == upload_option:

        #markdown css hack to remove fullscreen -- fickle because it is hardcoded
        nth_child_num += 1

        uploaded_file = st.file_uploader("Choose a file:", 
            type = rsc.valid_fasta_exts + rsc.valid_genbank_exts)

        if uploaded_file is not None:
            name, ext = rsc.get_name_ext(uploaded_file.name) # unused name -- could add in

            text_io = io.TextIOWrapper(uploaded_file, encoding='UTF-8')
            text = text_io.read() #saves this from losing in memory when stream is read
            st.success("File uploaded.")

            inSeq = rsc.validate_file(io.StringIO(text), ext)

    elif option == enter_option:

        inSeq = st.text_area('Input sequence here:',max_chars = rsc.maxPlasSize)
        inSeq = inSeq.replace("\n","")
        # remove numbers?
        rsc.validate_sequence(inSeq)
        
        #creates a procedurally-gen name for file based on seq 
        name = str(abs(hash(inSeq)))[:6]

    elif option == example_option:

        fastas=[]
        examples_path = rsc.get_example_fastas()
        for infile_loc in glob.glob(os.path.join(examples_path, "*.fa")):
            fastas.append(infile_loc.split("/")[-1].split(".fa")[0])
        exampleFile = st.radio("Choose example file:", fastas)
        inSeq = str(list(SeqIO.parse(os.path.join(examples_path, f"{exampleFile}.fa"), "fasta"))[0].seq)
        name = exampleFile

    if inSeq:

        with st.spinner("Annotating..."):
            linear = st.checkbox("Linear plasmid annotation")
            detailed = st.checkbox("Detailed plasmid annotation")

            with open(rsc.get_resource("templates", "FAQ.html")) as fh:
                faq = fh.read()
            sidebar.markdown(faq + images + cite_fund, unsafe_allow_html = True)

            recordDf = annotate(inSeq, args.blast_db, linear, detailed)

            if recordDf.empty:
                st.error("No annotations found.")
            else:
                recordDf = details(recordDf)
                st.markdown("---")
                st.header('Results:')

                st.write("Hover mouse for info, click and drag to pan, scroll wheel to zoom")
                st.bokeh_chart(get_bokeh(recordDf, linear), use_container_width=False)
                if linear:
                    st.write("\*plasmid is displayed as cirucular, though pLannotate is treating this as a linear construct")

                #markdown hack to remove full screen icon from bokeh plot (which is distorted)
                hide_full_screen = f'''
                <style>
                .element-container:nth-child({nth_child_num}) button{{visibility: hidden;}}
                </style>
                '''
                st.markdown(hide_full_screen, unsafe_allow_html=True)

                st.header("Download Annotations:")

                #write and encode gbk for dl
                gbk = rsc.get_gbk(recordDf, inSeq, linear)
                b64 = base64.b64encode(gbk.encode()).decode()
                gbk_dl = f'<a href="data:text/plain;base64,{b64}" download="{name}_pLann.gbk"> download {name}_pLann.gbk</a>'
                st.markdown(gbk_dl, unsafe_allow_html=True)

                #encode csv for dl
                cleaned = rsc.get_clean_csv_df(recordDf)
                csv = cleaned.to_csv(index=False)
                b64 = base64.b64encode(csv.encode()).decode()
                csv_dl = f'<a href="data:text/plain;base64,{b64}" download="{name}_pLann.csv"> download {name}_pLann.csv</a>'
                st.markdown(csv_dl, unsafe_allow_html=True)

                if option == upload_option and ext in rsc.valid_genbank_exts:
                    st.header("Download Combined Annotations:")
                    st.subheader("uploaded Genbank + pLannotate")
                    submitted_gbk = list(SeqIO.parse(io.StringIO(text), "gb"))[0]
                    gbk = rsc.get_gbk(recordDf, inSeq, linear, submitted_gbk)
                    b64 = base64.b64encode(gbk.encode()).decode()
                    gbk_dl = f'<a href="data:text/plain;base64,{b64}" download="{name}_pLann.gbk"> download {name}_pLann.gbk</a>'
                    st.markdown(gbk_dl, unsafe_allow_html=True)

                st.markdown("---")

                #prints table of features
                st.header("Features")
                displayColumns = ['Feature','percent identity','percent match length','Description',"database"]
                markdown = cleaned[displayColumns].copy()
                numericCols = ['percent identity', 'percent match length']
                markdown[numericCols] = np.round(markdown[numericCols], 1)
                markdown[numericCols] = markdown[numericCols].astype(str) + "%"
                markdown.loc[markdown['database'] == "infernal", 'percent identity'] = "-" #removes percent from infernal hits 
                markdown.loc[markdown['database'] == "infernal", 'percent match length'] = "-" #removes percent from infernal hits
                markdown = markdown.set_index("Feature",drop=True)
                markdown = markdown.drop("database", axis=1)
                st.markdown(markdown.drop_duplicates().to_markdown())


def setup_page():
    st.set_page_config(page_title="pLannotate", page_icon=rsc.get_image("icon.png"), layout='centered', initial_sidebar_state='auto')
    sys.tracebacklimit = 10 #removes traceback so code is not shown during errors

    hide_streamlit_style = '''
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    </style>
    '''
    #markdown hack to remove full screen icon from pLannotate logo
    hide_full_screen = '''
    <style>
    .element-container:nth-child(2) button{visibility: hidden;}
    </style>
    '''
    st.markdown(hide_full_screen, unsafe_allow_html=True)
    st.markdown(hide_streamlit_style, unsafe_allow_html=True)

    st.image(rsc.get_image("pLannotate.png"), use_column_width = "auto")

    st.markdown(f'<div style="text-align: right; font-size: 0.9em"> {plannotate_version} </div>', unsafe_allow_html=True)

    st.subheader('Annotate your engineered plasmids')
    sidebar = st.sidebar.empty()

    with open(rsc.get_template("blurb.html")) as fh:
        blurb = fh.read()
    with open(rsc.get_template("citation_funding.html")) as fh:
        cite_fund = fh.read()

    # have to use this b64-encoding hack to dislpay
    # local images, because html pathing is wonky/
    # not possible on streamlit
    with open(rsc.get_image("twitter.b64"), "r") as fh:
        twitter = fh.read()
    with open(rsc.get_image("email.b64"), "r") as fh:
        email = fh.read()
    with open(rsc.get_image("github.b64"), "r") as fh:
        github = fh.read()
    with open(rsc.get_image("paper.b64"), "r") as fh:
        paper = fh.read()

    # this is a python f-string saved as a .txt file
    # when processed, it becomes functional HTML
    # more readable than explicitly typing the f-string here
    newline = "\n"
    with open(rsc.get_resource("templates","images.txt")) as file:
        images = f"{file.read().replace(newline, '')}".format(
            twitter = twitter, email = email, github = github, paper = paper)

    sidebar.markdown(blurb + images + cite_fund, unsafe_allow_html=True)
    
    return sidebar,cite_fund,images
