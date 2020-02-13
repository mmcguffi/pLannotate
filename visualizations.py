from dna_features_viewer import BiopythonTranslator

def plot_plas(inRecord):
    class MyCustomTranslator(BiopythonTranslator):

        def compute_feature_color(self, feature):
            if "ColE1" in feature.qualifiers['label']:
                return "#4e7fff" #blue
            elif feature.type == "source":
                return "#4e7fff" #blue
            elif "origin" in feature.type :
                return "#4e7fff" #blue

            elif "Fragment" in feature.type :
                return "#808080" #grey

            elif "AmpR" in feature.qualifiers['label'] or "CmR" in feature.qualifiers['label']:
                return "#f6a35e" #orange

            elif feature.type == "CDS":
                return "#479f71" #green

            else:
                return "white"
        pass

    graphic_record = MyCustomTranslator().translate_record(inRecord,"circular")
    ax, _ = graphic_record.plot(figure_width=5)
    ax.figure.tight_layout()
    return ax

    # from PIL import Image
    # #removes extra whitespace at top of image. annoying hack
    # tempPic=NamedTemporaryFile(suffix='.png')
    # ax.figure.savefig(tempPic.name, bbox_inches="tight",transparent = True,)
    # img = Image.open(tempPic.name)
    # width, height = img.size
    # cropped = img.crop((0, 185, width, height-25))
    # cropped.save(tempPic.name)
    # st.image(tempPic.name)
    # tempPic.close()
