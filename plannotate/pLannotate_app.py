import base64
import glob
import io
import sys
from tempfile import NamedTemporaryFile

from Bio import SeqIO
import numpy as np
import pandas as pd
import streamlit as st

from plannotate import __version__ as plannotate_version
from plannotate.annotate import annotate, get_gbk
from plannotate.bokeh_plot import get_bokeh
from plannotate.BLAST_hit_details import details

def get_resource(group, name):
    return pkg_resources.resource_filename(__package__, f"{group}/{name}")

def get_image(name):
    return get_resource("images", name)

def get_template(name):
    return get_resource("templates", name)

def get_example_fastas():
    return get_resource("fastas", "")



def main():
    parser = argparse.ArgumentParser('pLannotate', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--version', action='version',
        version='%(prog)s {}'.format(__version__))
    parser.add_argument(
        "--blast_db", default="./BLAST_dbs/", help="path to BLAST databases.")
    parser.parse_args()

    st.set_page_config(page_title="pLannotate", page_icon=get_image("icon.png"), layout='centered', initial_sidebar_state='auto')
    sys.tracebacklimit = 0 #removes traceback so code is not shown during errors

    hide_streamlit_style = """
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    </style>
    """
    st.markdown(hide_streamlit_style, unsafe_allow_html=True)

    st.image(get_image("pLannotate.png"), use_column_width=False, width=500)

    #markdown hack to remove full screen icon from pLannotate logo
    hide_full_screen = '''
    <style>
    .element-container:nth-child(2) button{visibility: hidden;}
    </style>
    '''
    st.markdown(hide_full_screen, unsafe_allow_html=True)

    st.markdown(f'<div style="text-align: right; font-size: 0.9em"> {plannotate_version} </div>', unsafe_allow_html=True)

    st.subheader('Annotate your engineered plasmids')
    sidebar = st.sidebar.empty()

    with open(get_template("blurb.html")) as fh:
        blurb = fh.read()
    with open(get_template("citation_funding.html")) as fh:
        cite_fund = fh.read()

    # have to use this b64-encoding hack to dislpay
    # local images, because html pathing is wonky/
    # not possible on streamlit
    with open(get_image("twitter.b64"), "r") as fh:
        twitter = fh.read()
    with open(get_image("email.b64"), "r") as fh:
        email = fh.read()
    with open(get_image("github.b64"), "r") as fh:
        github = fh.read()
    with open(get_image("paper.b64"), "r") as fh:
        paper = fh.read()

    images = f'''
    <style>
        #images {{
        position: relative;
        bottom: -10px;
        left: 15px;
        }}
    </style>

    <div id='images'>
        <a href="https://twitter.com/matt_mcguffie">
            <img src="{twitter}"/>
        </a>
        <a href="mailto: mmcguffie@utexas.edu">
            <img src="{email}"/>
        </a>
        <a href="https://github.com/barricklab/pLannotate">
            <img src="{github}"/>
        </a>
        <a href="https://doi.org/10.1093/nar/gkab374">
            <img src="{paper}"/>
        </a>
    </div>
    <br>
    '''

    sidebar.markdown(blurb + images + cite_fund, unsafe_allow_html=True

    inSeq=""
    maxPlasSize = 50000
    IUPAC= 'GATCRYWSMKHBVDNgatcrywsmkhbvdn'

    # markdown css hack to remove fullscreen
    # fickle because it is hardcoded and can
    # change with streamlit versions updates
    nth_child_num = 14

    option = st.radio(
        'Choose method of submitting sequence:',
        ["Upload a file (.fa .fasta .gb .gbk)",
            "Enter a sequence",
            "Example"]
        )

    if option == "Upload a file (.fa .fasta .gb .gbk)":

        #markdown css hack to remove fullscreen -- fickle because it is hardcoded
        nth_child_num += 1

        uploaded_file = st.file_uploader("Choose a file:", type=['fa',"fasta","gb","gbk","gbff"])

        if uploaded_file is not None:
            extention = uploaded_file.name.split(".")[-1]
            if extention == "fa":
                extention = "fasta"
            elif extention == "gb" or extention == "gbff":
                extention = "gbk"

            text_io = io.TextIOWrapper(uploaded_file,encoding='UTF-8')

            st.success("File uploaded.")

            if extention == "fasta":

                #This catches errors on file uploads via Biopython
                fileloc = NamedTemporaryFile()
                record = list(SeqIO.parse(text_io, "fasta"))
                try:
                    record[0].annotations["molecule_type"] = "DNA"
                except IndexError:
                    error = "Malformed fasta file --> please submit a fasta file in standard format"
                    raise ValueError(error)
                SeqIO.write(record, fileloc.name, 'fasta')
                record = list(SeqIO.parse(fileloc.name, "fasta"))
                fileloc.close()

                if len(record)!=1:
                    error = 'FASTA file contains many entries --> please submit a single FASTA file.'
                    raise ValueError(error)

                inSeq = str(record[0].seq)

            elif extention == "gbk":
                fileloc = NamedTemporaryFile()
                try:
                    record = list(SeqIO.parse(text_io, "gb"))[0]
                except IndexError:
                    error = "Malformed Genbank file --> please submit a Genbank file in standard format"
                    raise ValueError(error)
                submitted_gbk = record
                SeqIO.write(record, fileloc.name, 'fasta')
                record = list(SeqIO.parse(fileloc.name, "fasta"))
                fileloc.close()
                inSeq = str(record[0].seq)

            else:
                st.error("invalid submission")

    elif option == "Enter a sequence":

        inSeq = st.text_area('Input sequence here:',max_chars = maxPlasSize)
        inSeq = inSeq.replace("\n","")

    elif option == "Example":

        fastas=[]
        examples_path = get_example_fastas()
        for infile_loc in glob.glob(os.path.join(examples_path, "*.fa")):
            fastas.append(infile_loc.split("/")[-1].split(".fa")[0])
        exampleFile = st.radio("Choose example file:", fastas)
        inSeq = str(list(SeqIO.parse(os.path.join(examples_path, f"{exampleFile}.fa"), "fasta"))[0].seq)

    if inSeq:

        if not set(inSeq).issubset(IUPAC):
            error = f'Sequence contains invalid characters -- must be ATCG and/or valid IUPAC nucleotide ambiguity code'
            raise ValueError(error)

        if len(inSeq) > maxPlasSize:
            error = f'Are you sure this is an engineered plasmid? Entry size is too large -- must be {maxPlasSize} bases or less.'
            raise ValueError(error)

        with st.spinner("Annotating..."):
            linear = st.checkbox("Linear plasmid annotation")

            with open("./FAQ.html") as fh:
                faq = fh.read()
            sidebar.markdown(faq + images + cite_fund, unsafe_allow_html = True)

            recordDf = annotate(inSeq, args.blast_db, linear)

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

                #creates a procedurally-gen name for file based on seq 
                filename = str(abs(hash(inSeq)))[:6]

                #write and encode gbk for dl
                gbk=get_gbk(recordDf, inSeq, linear)
                b64 = base64.b64encode(gbk.encode()).decode()
                gbk_dl = f'<a href="data:text/plain;base64,{b64}" download="{filename}.gbk"> download {filename}.gbk</a>'
                st.markdown(gbk_dl, unsafe_allow_html=True)

                #encode csv for dl
                columns = ['qstart', 'qend', 'sframe', 'pident', 'slen', 'sseq', 'length', 'uniprot', 'abs percmatch', 'fragment', 'db', 'Feature', 'Type', 'Description']
                cleaned = recordDf[columns]
                replacements = {'qstart':'start location', 'qend':'end location', 'sframe':'strand', 'pident':'percent identity', 'slen':'full length of feature in db', 'sseq':'full sequence of feature in db', 'length':'length of found feature', 'uniprot':'uniprot ID', 'abs percmatch':'percent match length', 'db':'database'}
                cleaned = cleaned.rename(columns=replacements)
                csv = cleaned.to_csv(index=False)
                b64 = base64.b64encode(csv.encode()).decode()
                csv_dl = f'<a href="data:text/plain;base64,{b64}" download="{filename}.csv"> download {filename}.csv</a>'
                st.markdown(csv_dl, unsafe_allow_html=True)

                if option == "Upload a file (.fa .fasta .gb .gbk)" and extention == "gbk":
                    st.header("Download Combined Annotations:")
                    st.subheader("uploaded Genbank + pLannotate")
                    gbk = get_gbk(recordDf, inSeq, linear, submitted_gbk)
                    b64 = base64.b64encode(gbk.encode()).decode()
                    gbk_dl = f'<a href="data:text/plain;base64,{b64}" download="{filename}.gbk"> download {filename}.gbk</a>'
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
