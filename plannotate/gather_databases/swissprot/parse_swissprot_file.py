# /bin/python3
import argparse
import re
from io import StringIO

import pandas as pd
from Bio import SeqIO


def swissprot(file_name, verbose):
    file = f"./split/{file_name}"

    # open file
    with open(file, "r") as handle:
        line = (
            handle.read() + "//"
        )  # because I trimmed this, though it is part of the specs
    with StringIO(line) as fastq_io:
        record_dict = SeqIO.to_dict(SeqIO.parse(fastq_io, "swiss"))

    swiss = record_dict[list(record_dict.keys())[0]]

    # this is all wonky and probably over-engineered
    id = swiss.id
    attrs = [
        "annotations,gene_name",
        "description",
        "annotations,comment",
        "annotations,protein_existence",
    ]
    details = {"name": swiss.name}
    for attr in attrs:
        attr = attr.split(",")
        try:
            if len(attr) == 1:
                details[attr[1]] = getattr(swiss, attr[0])
            elif len(attr) == 2:
                details[attr[1]] = getattr(swiss, attr[0])[attr[1]]
        except:
            pass

    try:
        protein_existence = pd.read_csv(
            "./protein_existence.csv", header=None, index_col=0, names=["description"]
        )
        protein_existence_level = details["protein_existence"]
        protein_existence_d = (
            f"{protein_existence.loc[protein_existence_level]['description']}: "
            f"Swiss-Prot protein existence level {protein_existence_level}. "
        )
    except (KeyError, IndexError):
        protein_existence = ""
        protein_existence_d = ""

    # pd.DataFrame([ele.split(":",1) for ele in details['comment'].split("\n")]).set_index(0).T
    try:
        comment = [ele.split(":", 1) for ele in details["comment"].split("\n")]
        comment = {ele[0]: ele[1].strip() for ele in comment}
    except KeyError:
        comment = {}

    try:
        function = comment["FUNCTION"]
        function_d = f"{function} "  # gives a space for stuff after
    except KeyError:
        function = ""
        function_d = ""

    try:
        biotech = comment["BIOTECHNOLOGY"]
        biotech_d = str(comment["BIOTECHNOLOGY"])  # copies the string
        # cuts out super long details (eg GFP)
        if len(biotech) > 200:
            biotech = ""
            biotech_d = ""
    except KeyError:
        biotech = ""
        biotech_d = ""

    try:
        organism = swiss.annotations["organism"]
        organism_d = f"From {organism}. "
    except KeyError:
        organism = ""
        organism_d = ""

    # tries to get a common name and an alt name, if available
    try:
        names = details["gene_name"].split(";")

        name = [ele for ele in names if "Name=" in ele][0]
        name = name.replace("Name=", "")

        swissprotName = swiss.name
        swissprotName_d = f"{swissprotName} - "  # gives a space for stuff after

    except (KeyError, IndexError):
        # altName = None
        name = details["name"]
        swissprotName = ""
        swissprotName_d = ""

    try:
        names = details["gene_name"].split(";")
        altName = [ele for ele in names if "Synonyms=" in ele][0]
        altName = altName.replace("Synonyms=", "")
        altName_d = f"Also known as {altName}. "
    except (KeyError, IndexError):
        altName = ""
        altName_d = ""

    anno = f"{swissprotName_d}{altName_d}{protein_existence_d}{function_d}{organism_d}{biotech_d}".replace(
        "  ", " "
    )
    anno = re.sub(r"\s*{.*}\s*.", "", anno).strip()  # removes text between curly braces
    name = re.sub(r"\s*{.*}\s*", " ", name).strip()  # this removes pubmed citations and
    # extra details on some names
    print(swissprotName)
    swissprotName = re.sub(r"\s*{.*}\s*.", "", swissprotName).strip()
    function = re.sub(r"\s*{.*}\s*.", "", function).strip()

    if verbose == False:
        swiss_description = pd.DataFrame(
            [
                id,
                name,
                swissprotName,
                altName,
                protein_existence_level,
                function,
                organism,
                biotech,
            ]
        ).T
    elif verbose == True:
        swiss_description = pd.DataFrame([id, name, anno]).T
    else:
        raise ValueError("verbose must be True or False")

    swiss_description.to_csv(
        "./swiss_description.csv", mode="a", header=False, index=False
    )

    return


# comand line arguments
parser = argparse.ArgumentParser(
    description="This script parses a swissprot file and returns a list of the proteins in the file"
)
parser.add_argument("file", help="The name of the file to parse")
# add a flag for verbose output
parser.add_argument("-v", "--verbose", help="fully formatted text", action="store_true")


args = parser.parse_args()
swissprot(args.file, args.verbose)
