"""Interactive Streamlit web app for pLannotate.

This is an optional front end (install with ``pip install 'plannotate[server]'``)
launched by ``plannotate streamlit``. It renders the same interface as before, but
builds its annotations on :class:`plannotate.models.Construct` rather than the legacy
``resources`` helpers, and embeds the Bokeh map as HTML so it works with Bokeh 3.

Streamlit runs this file as a script (``__name__ == "__main__"``), so it uses absolute
imports and only renders the page under that guard.
"""

# Streamlit is an optional ('server' extra) dependency; type checkers run in
# environments without it, so do not flag the import as unresolved.
# pyright: reportMissingImports=false, reportMissingModuleSource=false
from __future__ import annotations

import base64
import glob
import io
import os
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import streamlit as st
import streamlit.components.v1 as components
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from streamlit.delta_generator import DeltaGenerator

from plannotate import __version__, _package_data, validation
from plannotate.models import Construct

UPLOAD_OPTION = "Upload a file (FASTA or GenBank)"
ENTER_OPTION = "Enter a sequence"
EXAMPLE_OPTION = "Example"

FEATURE_TABLE_COLUMNS = [
    "Feature",
    "percent identity",
    "percent match length",
    "Description",
    "database",
]


def _yaml_file() -> Path:
    """Resolve the active database configuration set by ``plannotate streamlit``."""
    configured = os.environ.get("PLANNOTATE_YAML_FILE")
    return Path(configured) if configured else _package_data.get_yaml_path()


def _download_link(content: str, file_name: str) -> str:
    """Build a base64 data-URI download link, matching the original UI."""
    encoded = base64.b64encode(content.encode()).decode()
    return (
        f'<a href="data:text/plain;base64,{encoded}" '
        f'download="{file_name}"> download {file_name}</a>'
    )


def _setup_page() -> tuple[DeltaGenerator, str, str]:
    """Render the static page chrome and return the sidebar and its content blocks."""
    st.set_page_config(
        page_title="pLannotate",
        page_icon=str(_package_data.get_image("icon.png")),
        layout="centered",
        initial_sidebar_state="auto",
    )
    sys.tracebacklimit = 0  # hide tracebacks so source is not shown on errors

    st.markdown(
        "<style>#MainMenu {visibility: hidden;} footer {visibility: hidden;}</style>",
        unsafe_allow_html=True,
    )
    st.image(str(_package_data.get_image("pLannotate.png")), use_container_width=True)
    st.markdown(
        f'<div style="text-align: right; font-size: 0.9em"> {__version__} </div>',
        unsafe_allow_html=True,
    )
    st.subheader("Annotate your engineered plasmids")
    sidebar = st.sidebar.empty()

    blurb = _package_data.get_template("blurb.html").read_text()
    cite_fund = _package_data.get_template("citation_funding.html").read_text()

    # local images are inlined as base64 because Streamlit cannot serve file paths
    icons = {
        name: _package_data.get_image(f"{name}.b64").read_text()
        for name in ("twitter", "email", "github", "paper")
    }
    images_template = _package_data.get_template("images.txt").read_text()
    images = images_template.replace("\n", "").format(**icons)

    sidebar.markdown(blurb + images + cite_fund, unsafe_allow_html=True)

    st.markdown(
        "<style>"
        'button[title="View fullscreen"], button[title="View fullscreen"]:hover '
        "{display: none;}"
        "</style>",
        unsafe_allow_html=True,
    )
    return sidebar, cite_fund, images


def _read_record(text: str, ext: str) -> SeqRecord:
    """Parse and validate a single uploaded FASTA or GenBank record."""
    fmt = "genbank" if ext in validation.VALID_GENBANK_EXTS else "fasta"
    records = list(SeqIO.parse(io.StringIO(text), fmt))
    if len(records) != 1:
        raise validation.InvalidSequenceError(
            "File contains multiple entries -- please submit a single sequence file."
        )
    validation.validate_sequence(str(records[0].seq), max_length=None)
    return records[0]


def _collect_input() -> tuple[str, str, SeqRecord | None]:
    """Return (sequence, name, prior_record) for the chosen input method.

    ``prior_record`` is the uploaded GenBank record whose original features should be
    combined with pLannotate's, or None for FASTA / pasted / example input.
    """
    option = st.radio(
        "Choose method of submitting sequence:",
        [UPLOAD_OPTION, ENTER_OPTION, EXAMPLE_OPTION],
    )

    if option == UPLOAD_OPTION:
        uploaded = st.file_uploader(
            "Choose a file:",
            type=[
                ext.lstrip(".")
                for ext in validation.VALID_FASTA_EXTS + validation.VALID_GENBANK_EXTS
            ],
        )
        if uploaded is None:
            return "", "", None
        text = io.TextIOWrapper(uploaded, encoding="UTF-8").read()
        st.success("File uploaded.")
        name, ext = validation.get_name_ext(uploaded.name)
        record = _read_record(text, ext)
        prior = record if ext in validation.VALID_GENBANK_EXTS else None
        return str(record.seq), name, prior

    if option == ENTER_OPTION:
        entered = st.text_area(
            "Input sequence here:", max_chars=validation.MAX_PLAS_SIZE
        )
        sequence = "".join(char for char in entered if not char.isspace())
        sequence = "".join(char for char in sequence if not char.isdigit())
        if not sequence:
            return "", "", None
        validation.validate_sequence(sequence, max_length=None)
        return sequence, str(abs(hash(sequence)))[:6], None

    examples_path = _package_data.get_example_fastas()
    names = sorted(
        Path(path).stem for path in glob.glob(os.path.join(str(examples_path), "*.fa"))
    )
    chosen = st.radio("Choose example file:", names)
    record = SeqIO.read(os.path.join(str(examples_path), f"{chosen}.fa"), "fasta")
    return str(record.seq), chosen, None


def _feature_table(construct: Construct) -> str:
    """Render the feature summary table exactly as the original app did."""
    cleaned = construct.to_csv()
    table = cleaned[FEATURE_TABLE_COLUMNS].copy()
    percent_columns = ["percent identity", "percent match length"]
    table[percent_columns] = np.round(table[percent_columns], 1).astype(str) + "%"
    # Rfam hits have no meaningful percent identity / match length
    table.loc[table["database"] == "Rfam", percent_columns] = "-"
    table = table.set_index("Feature", drop=True).drop("database", axis=1)
    return table.drop_duplicates().to_markdown()


def _genbank(construct: Construct, record: SeqRecord | None) -> str:
    """Export GenBank, optionally layering pLannotate's features onto a base record."""
    return construct.to_genbank(record)


def _render_results(
    construct: Construct,
    name: str,
    linear: bool,
    detailed: bool,
    prior: SeqRecord | None,
) -> None:
    """Show the plasmid map, download links, and feature table for a construct."""
    st.markdown("---")
    st.header("Results:")
    st.write("Hover mouse for info, click and drag to pan, scroll wheel to zoom")
    components.html(construct.to_html(), height=820, width=820)

    if linear:
        st.write(
            r"\*plasmid is displayed as circular, though pLannotate is treating "
            r"this as a linear construct"
        )
    if detailed:
        st.write(
            r"\*\*pLannotate is running in Detailed Annotation mode which can find "
            r"more hits, though may also find more false positives."
        )

    st.header("Download Annotations:")
    # for a GenBank upload, keep its original header but only pLannotate's features
    annotations_only = prior
    if prior is not None:
        annotations_only = deepcopy(prior)
        annotations_only.features = []
    st.markdown(
        _download_link(_genbank(construct, annotations_only), f"{name}_pLann.gbk"),
        unsafe_allow_html=True,
    )
    csv = construct.to_csv().to_csv(index=False)
    st.markdown(
        _download_link(csv, f"{name}_pLann.csv"),
        unsafe_allow_html=True,
    )

    if prior is not None:
        st.header("Download Combined Annotations:")
        st.subheader("uploaded Genbank + pLannotate")
        st.markdown(
            _download_link(_genbank(construct, prior), f"{name}_pLann.gbk"),
            unsafe_allow_html=True,
        )

    st.markdown("---")
    st.header("Features")
    st.markdown(_feature_table(construct))


def render() -> None:
    """Render the full pLannotate web page."""
    sidebar, cite_fund, images = _setup_page()
    sequence, name, prior = _collect_input()
    if not sequence:
        return

    linear = st.checkbox("Linear plasmid annotation")
    detailed = st.checkbox("Detailed plasmid annotation")

    faq = _package_data.get_template("FAQ.html").read_text()
    sidebar.markdown(faq + images + cite_fund, unsafe_allow_html=True)

    with st.spinner("Annotating..."):
        construct = Construct(
            seq=sequence,
            linear=linear,
            detailed=detailed,
            db_options=_yaml_file(),
            name=name,
        )

    if not construct.features:
        st.error("No annotations found.")
        return

    _render_results(construct, name, linear, detailed, prior)


if __name__ == "__main__":
    render()
