"""Tests for the curated selection-marker and copy-number lookup tables."""

import pandas as pd
import pytest
import yaml

from plannotate import _curation, _package_data

MARKER_CLASSES = {
    "antibiotic resistance",
    "drug resistance",
    "herbicide resistance",
    "auxotrophic",
    "metabolic",
    "counter-selection",
    "screening",
}
COPY_CLASSES = {
    "very low",
    "low",
    "medium",
    "medium to high",
    "high",
    "not applicable",
    # no published figure supports even a class for this origin
    "unreported",
}
DOMAINS = {"bacterial", "eukaryotic", "both"}


def test_selection_marker_lookup():
    marker = _curation.selection_marker("snapgene", "AmpR")

    assert marker is not None
    assert marker.marker_class == "antibiotic resistance"
    assert "ampicillin" in marker.selection_agent
    assert marker.domain == "bacterial"
    assert marker.host_range == "broad (Gram-negative bacteria)"
    assert marker.reference == "PMID 358200"


def test_origin_copy_number_lookup():
    origin = _curation.origin_copy_number("snapgene", "p15A_ori")

    assert origin is not None
    assert origin.copy_number == "~20"
    assert origin.copy_class == "medium"
    assert "PMID" in origin.reference


def test_uncurated_records_return_nothing():
    assert _curation.selection_marker("snapgene", "EGFP") is None
    assert _curation.origin_copy_number("snapgene", "EGFP") is None
    # the two tables are disjoint: a marker is not an origin and vice versa
    assert _curation.origin_copy_number("snapgene", "AmpR") is None
    assert _curation.selection_marker("snapgene", "ori") is None


def test_global_feature_suppressions_are_source_pinned():
    suppressed = _curation.suppressed_feature_accessions()

    assert ("swissprot", "P03851") in suppressed
    assert ("snapgene", "ISS") in suppressed
    assert ("snapgene", "P03851") not in suppressed
    assert ("swissprot", "ISS") not in suppressed


def test_composite_reference_regions_are_source_pinned_and_coordinate_bearing():
    regions = _curation.composite_reference_regions()

    [t5_laco] = regions[("snapgene", "T5_promoter")]
    assert (t5_laco.start, t5_laco.end) == (21, 37)
    assert (t5_laco.component_db, t5_laco.component_sseqid) == (
        "snapgene",
        "lac_operator_(2)",
    )
    assert t5_laco.source == "curated_region:t5_laco_component"
    assert ("snapgene", "t5_promoter") not in regions


def test_fragment_suppression_regions_are_narrow_and_source_pinned():
    regions = _curation.fragment_suppression_regions()

    [pena] = regions[("swissprot", "Q02940")]
    assert (pena.subject_length, pena.start, pena.end) == (939, 580, 675)
    assert pena.max_identity == 70.0
    assert pena.source == "curated_fragment:puc_lac_region_penA"
    assert ("snapgene", "Q02940") not in regions


@pytest.mark.parametrize(
    ("updates", "duplicate", "message"),
    [
        ({"region_start": "not-an-integer"}, False, "must be integers"),
        ({"region_start": "5", "region_end": "4"}, False, "Invalid.*interval"),
        ({"source": ""}, False, "require pinned ids and provenance"),
        ({}, True, "Duplicate composite reference region"),
    ],
)
def test_malformed_composite_reference_regions_fail_loudly(
    tmp_path, monkeypatch, updates, duplicate, message
):
    row = {
        "db": "snapgene",
        "sseqid": "composite",
        "region_start": "1",
        "region_end": "3",
        "component_db": "snapgene",
        "component_sseqid": "component",
        "rationale": "test rationale",
        "source": "test:composite",
    } | updates
    path = tmp_path / "composite_reference_regions.csv"
    pd.DataFrame([row, row] if duplicate else [row]).to_csv(path, index=False)
    original = _package_data.get_resource

    def resource(group, filename):
        if filename == "composite_reference_regions.csv":
            return path
        return original(group, filename)

    monkeypatch.setattr(_package_data, "get_resource", resource)
    _curation.composite_reference_regions.cache_clear()
    try:
        with pytest.raises(ValueError, match=message):
            _curation.composite_reference_regions()
    finally:
        _curation.composite_reference_regions.cache_clear()


def test_an_accession_curated_for_one_database_does_not_match_another():
    # SnapGene keys its records by a name-derived slug and Swiss-Prot by accession, so
    # the same string can exist in one source and mean nothing in the other. SnapGene's
    # AMA1 is an Aspergillus replication origin; Swiss-Prot has three unrelated AMA1
    # records (an alpha-amanitin toxin, a Toxoplasma antigen, a yeast meiosis
    # regulator) that must not inherit the origin's claim.
    assert _curation.origin_copy_number("snapgene", "AMA1") is not None
    assert _curation.origin_copy_number("swissprot", "AMA1") is None
    assert _curation.selection_marker("snapgene", "DHFR") is not None
    assert _curation.selection_marker("swissprot", "P00374") is None
    # and a Swiss-Prot accession does not leak the other way
    assert _curation.selection_marker("swissprot", "P62593") is not None
    assert _curation.selection_marker("snapgene", "P62593") is None


def test_a_claim_holding_in_two_sources_is_two_rows():
    # ccdB is genuinely the same toxin in both sources, but an accession belongs to
    # one source only, so the claim is carried by one row per database.
    assert _curation.selection_marker("snapgene", "ccdB") is not None
    assert _curation.selection_marker("swissprot", "P62554") is not None


def test_pinning_splits_records_that_share_one_ambiguous_name():
    # Swiss-Prot has two records named `ura3`, and only one is the marker: P21594 is
    # PYRF (orotidine-5'-phosphate decarboxylase, the 5-FOA target) while P32747 is
    # PYRD (dihydroorotate dehydrogenase), which 5-FOA does not select against. Only
    # an accession key can tell them apart -- no name-level row could be right.
    good = _curation.selection_marker("swissprot", "P21594")
    assert good is not None
    assert "5-FOA" in good.selection_agent
    assert _curation.selection_marker("swissprot", "P32747") is None


def test_pinned_rows_exclude_the_entries_they_were_written_to_exclude():
    # `cat` names 17 chloramphenicol acetyltransferases and 3 CATA_* catalases
    assert _curation.selection_marker("swissprot", "P62577") is not None
    assert _curation.selection_marker("swissprot", "Q9PT92") is None
    # `sacB` names 5 levansucrases, a PTS protein, and 5 Neisseria capsule proteins
    assert _curation.selection_marker("swissprot", "P05655") is not None
    assert _curation.selection_marker("swissprot", "Q04938") is None
    assert _curation.selection_marker("swissprot", "Q9JWW8") is None


def test_known_trap_rows_pin_exactly_the_adjudicated_accessions():
    markers = pd.read_csv(
        _package_data.get_resource("data", "selection_markers.csv"), dtype=str
    ).set_index("name")

    assert markers.loc["cat", "sseqid"] == (
        "P00485;P00486;P06135;P23364;P50869;P26841;P36882;P62578;P62579;"
        "P62577;P58777;P07641;P62580;P25309;P20074;Q03058;P49417"
    )
    assert markers.loc["sacB", "sseqid"] == "P94468;P21130;P05655;F8DT26;P0DJA3"
    assert markers.loc["ura3", "sseqid"] == "P21594"
    assert "leu2" not in markers.index


def test_other_mixed_name_sets_exclude_unrelated_entries():
    # aadA Q4VR96 does not confer aminoglycoside resistance, ccdB Q2GPT7 is a
    # fungal biosynthetic enzyme, ccdB P45709 has no blurb support for the curated
    # gyrase-poison claim, and codA Q7X2H8 is choline oxidase.
    assert _curation.selection_marker("swissprot", "P0AG05") is not None
    assert _curation.selection_marker("swissprot", "Q4VR96") is None
    assert _curation.selection_marker("swissprot", "Q2GPT7") is None
    assert _curation.selection_marker("swissprot", "P45709") is None
    assert _curation.selection_marker("swissprot", "P25524") is not None
    assert _curation.selection_marker("swissprot", "Q7X2H8") is None


def test_lookup_keys_are_case_sensitive():
    assert _curation.selection_marker("snapgene", "AmpR") is not None
    assert _curation.selection_marker("snapgene", "ampr") is None
    assert _curation.selection_marker("SNAPGENE", "AmpR") is None


def test_lookup_tolerates_surrounding_whitespace():
    padded = _curation.selection_marker(" snapgene ", " KanR ")

    assert padded is not None
    assert padded == _curation.selection_marker("snapgene", "KanR")


@pytest.mark.parametrize(
    ("filename", "expected_columns"),
    [
        (
            "selection_markers.csv",
            [
                "db",
                "sseqid",
                "name",
                "marker_class",
                "selection_agent",
                "domain",
                "host_range",
                "reference",
            ],
        ),
        (
            "ori_copy_number.csv",
            [
                "db",
                "sseqid",
                "name",
                "copy_number",
                "copy_class",
                "domain",
                "host_range",
                "note",
                "reference",
            ],
        ),
        (
            "feature_suppressions.csv",
            ["db", "sseqid", "name", "rationale", "reference"],
        ),
        (
            "composite_reference_regions.csv",
            [
                "db",
                "sseqid",
                "region_start",
                "region_end",
                "component_db",
                "component_sseqid",
                "rationale",
                "source",
            ],
        ),
        (
            "fragment_suppression_regions.csv",
            [
                "db",
                "sseqid",
                "name",
                "subject_length",
                "region_start",
                "region_end",
                "max_identity",
                "rationale",
                "source",
            ],
        ),
    ],
)
def test_curated_tables_are_well_formed(filename, expected_columns):
    table = pd.read_csv(_package_data.get_resource("data", filename), dtype=str).fillna(
        ""
    )

    assert list(table.columns) == expected_columns
    # every row declares a source vocabulary; without one it could never be matched
    assert (table["db"].str.strip() != "").all()
    if "name" in table:
        assert (table["name"].str.strip() != "").all()
    for column in set(expected_columns) & {
        "db",
        "name",
        "sseqid",
        "component_db",
        "component_sseqid",
    }:
        assert (table[column] == table[column].str.strip()).all(), column

    # Every claim is frozen to the packaged source records that were manually checked.
    assert (table["sseqid"] != "").all(), f"{filename} has an unpinned row"

    # An accession belongs to exactly one source, so a row names exactly one database.
    # Pairing a ";"-separated db list with pinned accessions would take their cross
    # product and mint keys like (snapgene, <a Swiss-Prot accession>) that can never
    # match -- silently, since a key that matches nothing looks just like a rare hit.
    assert not table["db"].str.contains(";").any(), sorted(
        table.loc[table["db"].str.contains(";"), "name"]
    )

    # The accession join key must be unique per source.
    keys = [
        (row["db"].strip(), accession.strip())
        for _, row in table.iterrows()
        for accession in row["sseqid"].split(";")
    ]
    duplicates = {key for key in keys if keys.count(key) > 1}
    assert not duplicates, sorted(duplicates)


def test_curated_tables_use_the_documented_vocabularies():
    markers = pd.read_csv(
        _package_data.get_resource("data", "selection_markers.csv"), dtype=str
    )
    origins = pd.read_csv(
        _package_data.get_resource("data", "ori_copy_number.csv"), dtype=str
    )

    assert set(markers["marker_class"]) <= MARKER_CLASSES
    assert set(origins["copy_class"]) <= COPY_CLASSES
    # an origin with no published figure still gets a class, never a bare blank row
    assert origins["copy_class"].notna().all()

    for table in (markers, origins):
        assert set(table["domain"]) <= DOMAINS, sorted(set(table["domain"]) - DOMAINS)
        # host_range states its breadth first so the two columns cannot drift apart
        assert table["host_range"].str.startswith(("broad", "narrow")).all()

    # a eukaryote-only marker cannot be selected in a bacterial host, and vice versa
    bacterial = markers.loc[markers["domain"] == "bacterial", "host_range"]
    assert not bacterial.str.contains("mammalian|plant|yeast|fungi").any()


def test_curated_db_values_name_real_annotation_sources():
    # a typo in `db` is silent: the row simply never matches anything. The packaged
    # config is the authority on what a source can be called, so check against it.
    with _package_data.get_yaml_path().open() as handle:
        sources = set(yaml.safe_load(handle))

    used = {
        str(database).strip()
        for filename in (
            "selection_markers.csv",
            "ori_copy_number.csv",
            "feature_suppressions.csv",
            "composite_reference_regions.csv",
            "fragment_suppression_regions.csv",
        )
        for database in pd.read_csv(
            _package_data.get_resource("data", filename), dtype=str
        )["db"]
    }
    composite = pd.read_csv(
        _package_data.get_resource("data", "composite_reference_regions.csv"),
        dtype=str,
    )
    used.update(composite["component_db"].str.strip())

    assert used <= sources, sorted(used - sources)


def test_every_stated_copy_number_cites_a_reference():
    # the figures are strain- and growth-rate-dependent, so an uncited number is a
    # number nobody can check; an origin with no source states a class and nothing more
    origins = pd.read_csv(
        _package_data.get_resource("data", "ori_copy_number.csv"), dtype=str
    ).fillna("")

    uncited = origins[(origins["copy_number"] != "") & (origins["reference"] == "")]

    assert set(uncited["name"]) <= {"oriV", "ori2", "ARS1", "ars1"}, sorted(
        uncited["name"]
    )
