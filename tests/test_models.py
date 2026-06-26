"""Tests for public construct and feature models."""

from io import StringIO
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from plannotate import Construct as PublicConstruct
from plannotate import Feature as PublicFeature
from plannotate.models import Construct, Feature, df_to_features
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


def test_construct_plot_and_html_resources(annotated_construct):
    plot = annotated_construct.plot()
    cdn_html = annotated_construct.to_html()
    inline_html = annotated_construct.to_html(htmlfull=True)

    assert plot is not None
    assert "<!DOCTYPE html>" in cdn_html
    assert len(inline_html) > len(cdn_html)
