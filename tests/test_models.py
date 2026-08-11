"""Tests for public construct and feature models."""

from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from plannotate import Construct as PublicConstruct
from plannotate import Feature as PublicFeature
from plannotate.models import Construct, Feature, df_to_features, record_locus_name
from plannotate.validation import InvalidSequenceError

TEST_DATA = Path(__file__).parent / "test_data"


@pytest.fixture
def annotated_construct():
    annotations = pd.read_csv(TEST_DATA / "pXampl3.csv")
    sequence = SeqIO.read(TEST_DATA / "pXampl3.fa", "fasta").seq
    construct = Construct(sequence, detailed=True, _skip_annotation=True)
    construct.features = df_to_features(annotations)
    return construct


def test_domain_models_are_exposed_at_package_root():
    assert PublicConstruct is Construct
    assert PublicFeature is Feature


def test_construct_preserves_prior_annotations_without_mutating_input():
    prior = SeqRecord(Seq("ACGT"), id="existing", name="existing")
    prior.annotations["comment"] = "Original annotation"
    prior.features.append(SeqFeature(FeatureLocation(0, 2), type="promoter"))

    construct = Construct(
        "ACGT",
        prior_annotations=prior,
        _skip_annotation=True,
    )
    output = construct.to_seqrecord()
    comment = str(output.annotations["comment"])

    assert len(output.features) == 1
    assert output.features[0].type == "promoter"
    assert "Original annotation" in comment
    assert "pLannotate" in comment
    assert prior.annotations["comment"] == "Original annotation"


def test_construct_rejects_prior_annotations_for_another_sequence():
    prior = SeqRecord(Seq("AAAA"))

    with pytest.raises(ValueError, match="does not match"):
        Construct("ACGT", prior_annotations=prior, _skip_annotation=True)


@pytest.mark.parametrize("sequence", ["", "ACXT"])
def test_construct_rejects_invalid_sequences(sequence):
    with pytest.raises(InvalidSequenceError):
        Construct(sequence, _skip_annotation=True)


def test_construct_rejects_invalid_core_count():
    with pytest.raises(ValueError, match="cores must be at least 1"):
        Construct("ACGT", cores=0, _skip_annotation=True)


def test_construct_uses_prior_record_name_by_default():
    prior = SeqRecord(Seq("ACGT"), id="existing", name="existing")

    construct = Construct("ACGT", prior_annotations=prior, _skip_annotation=True)

    assert construct.name == "existing"


def test_html_export_does_not_display_plot(monkeypatch):
    construct = Construct("ACGT", _skip_annotation=True)

    monkeypatch.setattr(
        construct,
        "plot",
        lambda *args, **kwargs: pytest.fail("to_html() displayed the plot"),
    )

    assert "pLannotate" in construct.to_html()


def test_dataframe_feature_conversion():
    dataframe = pd.read_csv(TEST_DATA / "pXampl3.csv")

    features = df_to_features(dataframe)

    assert len(features) == len(dataframe)
    assert any(feature.is_forward_strand for feature in features)
    assert any(feature.is_reverse_strand for feature in features)
    assert df_to_features(pd.DataFrame()) == []


def test_every_annotation_column_maps_to_a_feature_field():
    from dataclasses import fields

    from plannotate._schema import ANNOTATION_COLUMNS
    from plannotate.models import _field

    feature_fields = {item.name for item in fields(Feature)}
    for column in ANNOTATION_COLUMNS:
        assert _field(column) in feature_fields, column


def test_construct_exports(annotated_construct):
    record = annotated_construct.to_seqrecord()
    genbank = SeqIO.read(StringIO(annotated_construct.to_genbank()), "genbank")
    csv = annotated_construct.to_csv()

    assert len(record.features) == 21
    assert len(genbank.features) == 21
    assert record.annotations["topology"] == "circular"
    assert "pLannotate" in record.annotations["comment"]
    assert len(csv.columns) == 14
    assert "start location" in csv.columns


def _feature(**overrides):
    values: dict[str, Any] = dict(
        sseqid="feat",
        feature_name="AmpR",
        description="beta-lactamase",
        feature_type="CDS",
        database="snapgene",
        qstart=0,
        qend=100,
        qlen=1000,
        sstart=1,
        send=101,
        sframe=1,
        qseq="A" * 100,
        length=100,
        slen=101,
        pident=100.0,
        percmatch=100.0,
        abs_percmatch=100.0,
        pi_permatch=100.0,
        evalue=0.0,
        score=100.0,
        priority=1,
        kind=1,
        fragment=False,
        wiggle=15,
        wstart=15,
        wend=85,
    )
    values.update(overrides)
    return Feature(**values)


def test_seqfeature_reports_the_matched_subject_range_and_traceback():
    qualifiers = _feature(sstart=12, send=112, btop="50AG49").seqfeature.qualifiers

    assert qualifiers["subject_start"] == 12
    assert qualifiers["subject_end"] == 112
    assert qualifiers["btop"] == "50AG49"


def test_seqfeature_omits_an_empty_traceback():
    # Infernal reports covariance-model hits with no base-by-base traceback
    assert "btop" not in _feature(btop="").seqfeature.qualifiers


@pytest.mark.parametrize(
    ("btop", "expected"),
    [(None, ""), (float("nan"), ""), (300, "300")],
)
def test_feature_normalizes_a_non_string_traceback(btop, expected):
    # a gapless, fully identical alignment is a bare match run, so a DataFrame column
    # of such hits can arrive numeric rather than as text
    assert _feature(btop=btop).btop == expected


def test_feature_defaults_the_traceback_for_a_csv_without_one():
    # CSVs written before btop existed must still round-trip through Feature
    legacy = pd.read_csv(TEST_DATA / "pXampl3.csv").drop(
        columns="btop", errors="ignore"
    )

    features = df_to_features(legacy)

    assert len(features) == len(legacy)
    assert all(feature.btop == "" for feature in features)


def test_seqfeature_annotates_a_selection_marker():
    qualifiers = _feature(sseqid="KanR", feature_name="KanR").seqfeature.qualifiers

    assert qualifiers["selection_marker"] == "antibiotic resistance"
    assert "kanamycin" in qualifiers["selection_agent"]
    assert qualifiers["domain"] == "both"
    assert qualifiers["host_range"] == "broad (bacteria and eukaryotes)"
    assert qualifiers["reference"] == "PMID 6270337"
    assert "copy_number" not in qualifiers


def test_seqfeature_annotates_an_origin_copy_number():
    qualifiers = _feature(
        sseqid="pSC101_ori", feature_name="pSC101 ori", feature_type="rep_origin"
    ).seqfeature.qualifiers

    assert qualifiers["copy_number"] == "~5"
    assert qualifiers["copy_number_class"] == "low"
    assert "Rep101" in qualifiers["copy_number_note"]
    assert qualifiers["reference"] == "PMID 29371642"
    assert qualifiers["domain"] == "bacterial"
    assert qualifiers["host_range"] == "narrow (enterobacteria)"
    assert "selection_marker" not in qualifiers


def test_seqfeature_omits_an_unknown_copy_number_but_keeps_its_class():
    qualifiers = _feature(
        sseqid="f1_ori", feature_name="f1 ori", feature_type="rep_origin"
    ).seqfeature.qualifiers

    assert "copy_number" not in qualifiers
    assert qualifiers["copy_number_class"] == "not applicable"


def test_seqfeature_curation_is_scoped_to_the_source_database():
    # DHFR is a methotrexate selection marker in SnapGene, but Swiss-Prot's DHFR
    # entries are ordinary dihydrofolate reductases -- the same label, a different
    # claim. The curated qualifiers must follow the database the hit came from.
    snapgene = _feature(
        sseqid="DHFR", feature_name="DHFR", database="snapgene"
    ).seqfeature
    swissprot = _feature(
        sseqid="P00374", feature_name="DHFR", database="swissprot"
    ).seqfeature

    assert snapgene.qualifiers["selection_marker"] == "drug resistance"
    assert "selection_marker" not in swissprot.qualifiers


def test_seqfeature_leaves_an_uncurated_feature_alone():
    qualifiers = _feature(feature_name="EGFP").seqfeature.qualifiers

    curated = {
        "selection_marker",
        "selection_agent",
        "copy_number",
        "copy_number_class",
        "copy_number_note",
        "domain",
        "host_range",
        "reference",
    }
    assert curated.isdisjoint(qualifiers)


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("lambda att insert seq", "lambda_att_insert_seq"),
        ("  spaced  out  ", "spaced_out"),
        ("with\ttab", "with_tab"),
        ("   ", "construct"),
        (None, "construct"),
    ],
)
def test_genbank_locus_name_tolerates_whitespace(name, expected):
    # the web app names constructs after the uploaded file, whose stem may contain
    # spaces; Biopython rejects whitespace on the LOCUS line
    construct = Construct("ACGT", _skip_annotation=True, name=name)

    record = SeqIO.read(StringIO(construct.to_genbank()), "genbank")

    assert record.name == expected


def test_genbank_locus_name_tolerates_whitespace_in_prior_record():
    prior = SeqRecord(Seq("ACGT"), id="prior name", name="prior name")

    construct = Construct("ACGT", _skip_annotation=True, prior_annotations=prior)
    record = SeqIO.read(StringIO(construct.to_genbank()), "genbank")

    assert record.name == "prior_name"
    assert prior.name == "prior name"  # the caller's record is not mutated


def test_genbank_keeps_biopython_placeholder_name():
    # Biopython special-cases its own placeholder into a bare "." locus, so
    # sanitizing it would replace valid output with a bogus name
    construct = Construct(
        "ACGT", _skip_annotation=True, prior_annotations=SeqRecord(Seq("ACGT"))
    )

    assert "LOCUS       .  " in construct.to_genbank()


@pytest.mark.parametrize(
    ("record", "is_genbank", "expected"),
    [
        # GenBank names itself from its LOCUS line, never from its accession
        (
            SeqRecord(Seq("ACGT"), id="AB123456.7", name="FriendlyLocus"),
            True,
            "FriendlyLocus",
        ),
        (SeqRecord(Seq("ACGT"), id="plasmidA", name="plasmidA"), False, "plasmidA"),
        # a bare ">" header parses to an empty id, naming nothing
        (SeqRecord(Seq("ACGT"), id="", name=""), False, None),
        # Biopython's placeholders are not names either
        (SeqRecord(Seq("ACGT")), False, None),
        (SeqRecord(Seq("ACGT")), True, None),
    ],
)
def test_record_locus_name_uses_the_field_each_format_names_itself_by(
    record, is_genbank, expected
):
    assert record_locus_name(record, is_genbank) == expected


def test_construct_plot_and_html_resources(annotated_construct):
    plot = annotated_construct.plot()
    cdn_html = annotated_construct.to_html()
    inline_html = annotated_construct.to_html(htmlfull=True)

    assert plot is not None
    assert "<!DOCTYPE html>" in cdn_html
    assert len(inline_html) > len(cdn_html)
