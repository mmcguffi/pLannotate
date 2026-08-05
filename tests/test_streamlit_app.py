"""Tests for the optional Streamlit front end's input handling."""

import io

import pytest

streamlit_app = pytest.importorskip(
    "plannotate.streamlit_app",
    reason="the web app needs the 'server' extra",
)


class _Upload(io.BytesIO):
    """Stand-in for Streamlit's uploaded-file object, which carries a name."""

    def __init__(self, name: str, text: str):
        super().__init__(text.encode())
        self.name = name


def _answer_prompts(monkeypatch, *answers):
    """Drive the app's radio prompts with canned answers, in order."""
    remaining = iter(answers)
    monkeypatch.setattr(
        streamlit_app.st, "radio", lambda *args, **kwargs: next(remaining)
    )


def _ignore_status_messages(monkeypatch):
    monkeypatch.setattr(streamlit_app.st, "success", lambda *args, **kwargs: None)


def test_example_input_is_named_after_its_file_not_its_header(monkeypatch):
    _answer_prompts(monkeypatch, streamlit_app.EXAMPLE_OPTION, "pCMVR8.74")

    _sequence, file_name, locus_name, prior = streamlit_app._collect_input()

    # this example's FASTA header reads ">Addgene", which is no one's plasmid name
    assert file_name == "pCMVR8.74"
    assert locus_name == "pCMVR8.74"
    assert prior is None


def test_uploaded_fasta_is_named_after_its_record(monkeypatch):
    _answer_prompts(monkeypatch, streamlit_app.UPLOAD_OPTION)
    _ignore_status_messages(monkeypatch)
    monkeypatch.setattr(
        streamlit_app.st,
        "file_uploader",
        lambda *args, **kwargs: _Upload(
            "some file name.fa", ">plasmidA\nACGTACGTACGT\n"
        ),
    )

    _sequence, file_name, locus_name, prior = streamlit_app._collect_input()

    # the file names the download, the record names the locus -- as in the CLI
    assert file_name == "some file name"
    assert locus_name == "plasmidA"
    assert prior is None


@pytest.mark.parametrize(
    ("header", "expected_locus"),
    [(">plasmidA description here\n", "plasmidA"), (">\n", "")],
)
def test_uploaded_fasta_locus_comes_from_the_header_id(
    monkeypatch, header, expected_locus
):
    _answer_prompts(monkeypatch, streamlit_app.UPLOAD_OPTION)
    _ignore_status_messages(monkeypatch)
    monkeypatch.setattr(
        streamlit_app.st,
        "file_uploader",
        lambda *args, **kwargs: _Upload("upload.fa", f"{header}ACGTACGTACGT\n"),
    )

    *_, locus_name, _prior = streamlit_app._collect_input()

    # an empty locus name defers to Construct's fallback, exactly as the CLI does
    assert locus_name == expected_locus


def test_uploaded_genbank_lets_its_own_record_name_the_locus(monkeypatch):
    genbank = (
        "LOCUS       FriendlyLocus              12 bp    DNA     circular UNK\n"
        "ACCESSION   AB123456\n"
        "VERSION     AB123456.7\n"
        "FEATURES             Location/Qualifiers\n"
        "ORIGIN\n"
        "        1 acgtacgtacgt\n"
        "//\n"
    )
    _answer_prompts(monkeypatch, streamlit_app.UPLOAD_OPTION)
    _ignore_status_messages(monkeypatch)
    monkeypatch.setattr(
        streamlit_app.st,
        "file_uploader",
        lambda *args, **kwargs: _Upload("upload.gbk", genbank),
    )

    _sequence, _file_name, locus_name, prior = streamlit_app._collect_input()

    # the LOCUS line names it, not the AB123456.7 accession
    assert locus_name == "FriendlyLocus"
    assert prior is not None
    assert prior.name == "FriendlyLocus"
