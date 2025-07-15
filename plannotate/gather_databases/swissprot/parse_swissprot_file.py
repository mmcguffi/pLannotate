# /bin/python3
import argparse
import os
import re
from io import StringIO

import pandas as pd
from Bio import SeqIO


def swissprot(file_name, output_file):
    file = file_name

    # Read raw file content to handle multiple records
    with open(file, "r") as handle:
        content = handle.read()

    # Split content into complete records (ending with //)
    # Use regex to split on // at the beginning of lines
    records_raw = re.split(r"\n//", content)

    for record_raw in records_raw:
        record_raw = record_raw.strip()
        if not record_raw or not record_raw.startswith("ID"):
            continue  # Skip empty or incomplete records

        # Add back the // terminator for proper parsing
        record_raw += "\n//"

        try:
            # Parse the individual record
            with StringIO(record_raw) as record_handle:
                record_dict = SeqIO.to_dict(SeqIO.parse(record_handle, "swiss"))
                swiss = record_dict[list(record_dict.keys())[0]]
                # Extract PE line manually since BioPython doesn't parse it
                pe_line = None
                for line in record_raw.split("\n"):
                    if line.startswith("PE"):
                        pe_line = line
                        break
                process_single_record(swiss, output_file, pe_line)
        except Exception as e:
            print(f"Warning: Could not parse record in {file_name}: {e}")
            continue

    return


def process_single_record(swiss, output_file, pe_line=None):
    """Process a single Swiss-Prot record"""
    try:
        # this is all wonky and probably over-engineered
        id = swiss.id
        attrs = [
            "annotations,gene_name",
            "description",
            "annotations,comment",
        ]
        details = {"name": swiss.name}
        for attr in attrs:
            attr_parts = attr.split(",")
            try:
                if len(attr_parts) == 1:
                    details[attr_parts[0]] = getattr(swiss, attr_parts[0])
                elif len(attr_parts) == 2:
                    details[attr_parts[1]] = getattr(swiss, attr_parts[0])[
                        attr_parts[1]
                    ]
            except Exception:
                pass

        # Extract protein existence level from PE line manually
        try:
            if pe_line:
                # Extract PE level from line like "PE   4: Predicted;"
                pe_match = re.search(r"PE\s+(\d+):", pe_line)
                if pe_match:
                    protein_existence_level = int(pe_match.group(1))
                    # Get path relative to this script
                    script_dir = os.path.dirname(os.path.abspath(__file__))
                    protein_existence = pd.read_csv(
                        os.path.join(script_dir, "swissprot-existance-levels.csv"),
                        header=None,
                        index_col=0,
                        names=["description"],
                    )
                    protein_existence_d = (
                        f"{protein_existence.loc[protein_existence_level]['description']}: "
                        f"Swiss-Prot protein existence level {protein_existence_level}. "
                    )
                else:
                    protein_existence_level = ""
                    protein_existence_d = ""
            else:
                protein_existence_level = ""
                protein_existence_d = ""
        except (KeyError, IndexError):
            protein_existence = ""
            protein_existence_d = ""
            protein_existence_level = ""

        # pd.DataFrame([ele.split(":",1) for ele in details['comment'].split("\n")]).set_index(0).T
        try:
            comment_data = details["comment"]
            # Handle case where comment might be a list (from BioPython)
            if isinstance(comment_data, list):
                # Convert list elements to strings, handling nested structures
                comment_text = "\n".join(str(item) for item in comment_data)
            else:
                comment_text = str(comment_data)
            comment = [
                ele.split(":", 1) for ele in comment_text.split("\n") if ":" in ele
            ]
            comment = {ele[0]: ele[1].strip() for ele in comment if len(ele) == 2}
        except (KeyError, TypeError):
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
            gene_name_data = details["gene_name"]

            # Handle BioPython's structured gene_name format
            if isinstance(gene_name_data, list):
                # Look for primary Name first
                name = None
                for gene_entry in gene_name_data:
                    if isinstance(gene_entry, dict) and "Name" in gene_entry:
                        if isinstance(gene_entry["Name"], list):
                            name = gene_entry["Name"][0]
                        else:
                            name = gene_entry["Name"]
                        # Clean up any ECO annotations
                        name = re.sub(r"\s*{.*?}\s*", "", name).strip()
                        break

                # If no Name found, try ORFNames
                if not name:
                    for gene_entry in gene_name_data:
                        if isinstance(gene_entry, dict) and "ORFNames" in gene_entry:
                            orf_names = gene_entry["ORFNames"]
                            if isinstance(orf_names, list) and orf_names:
                                name = orf_names[0]
                                # Clean up any ECO annotations
                                name = re.sub(r"\s*{.*?}\s*", "", name).strip()
                                break

                # If still no name found, use the swiss.name
                if not name:
                    name = details["name"]
            else:
                # Handle old string format if it exists
                gene_name_text = str(gene_name_data)
                names = gene_name_text.split(";")
                name_parts = [ele for ele in names if "Name=" in ele]
                if name_parts:
                    name = name_parts[0].replace("Name=", "").strip()
                    name = re.sub(r"\s*{.*?}\s*", "", name).strip()
                else:
                    name = details["name"]

            # swissprotName should be the Swiss-Prot entry name, not the gene name
            swissprotName = (
                swiss.name
            )  # This is the Swiss-Prot entry name (e.g., 001R_FRG3G)
            swissprotName_d = f"{swissprotName} - "  # gives a space for stuff after

        except (KeyError, IndexError, TypeError):
            # altName = None
            name = details["name"]
            swissprotName = ""
            swissprotName_d = ""

        try:
            gene_name_data = details["gene_name"]
            altName = ""
            altName_d = ""

            # Handle BioPython's structured gene_name format
            if isinstance(gene_name_data, list):
                for gene_entry in gene_name_data:
                    if isinstance(gene_entry, dict) and "Synonyms" in gene_entry:
                        synonyms = gene_entry["Synonyms"]
                        if isinstance(synonyms, list) and synonyms:
                            altName = synonyms[0]  # Take the first synonym
                            # Clean up any ECO annotations
                            altName = re.sub(r"\s*{.*?}\s*", "", altName).strip()
                            altName_d = f"Also known as {altName}. "
                            break
            else:
                # Handle old string format if it exists
                gene_name_text = str(gene_name_data)
                names = gene_name_text.split(";")
                synonym_parts = [ele for ele in names if "Synonyms=" in ele]
                if synonym_parts:
                    altName = synonym_parts[0].replace("Synonyms=", "").strip()
                    altName = re.sub(r"\s*{.*?}\s*", "", altName).strip()
                    altName_d = f"Also known as {altName}. "
        except (KeyError, IndexError, TypeError):
            altName = ""
            altName_d = ""

        anno = f"{swissprotName_d}{altName_d}{protein_existence_d}{function_d}{organism_d}{biotech_d}".replace(
            "  ", " "
        )
        anno = re.sub(
            r"\s*{.*}\s*.", "", anno
        ).strip()  # removes text between curly braces
        name = re.sub(
            r"\s*{.*}\s*", " ", name
        ).strip()  # this removes pubmed citations and
        # extra details on some names
        swissprotName = re.sub(r"\s*{.*}\s*.", "", swissprotName).strip()
        function = re.sub(r"\s*{.*}\s*.", "", function).strip()

        # Standard 4-column format: sseqid, name, type, blurb
        swiss_description = pd.DataFrame([id, name, "CDS", anno]).T

        swiss_description.to_csv(
            output_file, mode="a", header=False, index=False, sep='\t'
        )
    except Exception as e:
        # Skip records that can't be processed
        print(
            f"Warning: Could not process record {getattr(swiss, 'id', 'unknown')}: {e}"
        )
        pass


# command line arguments
parser = argparse.ArgumentParser(
    description="This script parses a swissprot file and returns a list of the proteins in the file"
)
parser.add_argument("file", help="The name of the file to parse")
parser.add_argument("-o", "--output", help="Output file path", required=True)


args = parser.parse_args()
swissprot(args.file, args.output)
