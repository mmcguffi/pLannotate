from bokeh.model import collect_filtered_models
import pandas as pd
import random

def gbk_to_df(gbk):
    features = []
    cols = ["qstart","qend","qlen","sframe","Feature","Description","Type","level","pi_permatch","fragment"]
    for feat in gbk.features:
        Type  = feat.type
        
        if Type == 'source':
            continue
        
        start = int(feat.location.start)
        end = int(feat.location.end)
        qlen = len(gbk.seq)
        strand = feat.location.strand

        try:
            Feature = feat.qualifiers['label'][0]
        except KeyError:
            Feature = "other"
        
        try:
            Description = feat.qualifiers['note'][0]
        except KeyError:
            Description = ""
        
        level = random.randint(0, 4)
        pi_permatch = 100
        fragment = False

        features.append((start,end,qlen,strand,Feature,Description,Type,level,pi_permatch,fragment))
    features = pd.DataFrame(features,columns=cols)
    return(features)

# "Upload an annotated Genbank file to add annotations (.gb or .gbk)",

    # elif option == "Upload an annotated Genbank file to add annotations (.gb or .gbk)":
    # #markdown css hack to remove fullscreen -- fickle because it is hardcoded
    # nth_child_num = 13

    # uploaded_file = st.file_uploader("Choose a file:", type=['gb',"gbk"])

    # if uploaded_file is not None:
    #     text_io = io.TextIOWrapper(uploaded_file,encoding='UTF-8')

    #     st.success("File uploaded.")
        
    #     record = list(SeqIO.parse(text_io, "gb"))[0]
    #     record = gbk_to_df(record)
    #     st.write(record)
    #     st.bokeh_chart(get_bokeh(record), use_container_width=False)
