"""Tests for the database nested-feature audit."""

from pathlib import Path

import pandas as pd

from plannotate._nested import (
    suppress_curated_fragment_artifacts,
    suppress_nested_fragments,
    suppress_uninformative_composite_fragments,
)
from tools.nested_feature_audit import (
    _align_metadata_ids,
    _back_translate,
    nested_rows,
)
from tools.nested_feature_figures import build_records
from tools.nested_feature_rules import (
    classify_report,
    classify_row,
    curated_decisions,
)

ROOT = Path(__file__).resolve().parents[1]


def _annotations(rows):
    defaults = {
        "name": "feature",
        "type": "CDS",
        "sframe": 1,
        "qstart": 0,
        "qend": 100,
        "length": 100,
        "slen": 100,
        "sstart": 1,
        "send": 100,
        "pident": 100.0,
        "abs percmatch": 100.0,
        "evalue": 1e-20,
        "fragment": False,
    }
    return pd.DataFrame([defaults | row for row in rows])


def test_align_metadata_ids_accepts_unique_case_only_drift():
    metadata = {"E4orf6": {"name": "E4orf6", "type": "CDS"}}

    assert _align_metadata_ids("snapgene", {"E4ORF6": "ATG"}, metadata) == {
        "E4ORF6": metadata["E4orf6"]
    }


def test_back_translate_preserves_protein_and_masks_unknown_residues():
    assert _back_translate("M*X") == "ATGTAANNN"


def test_nested_rows_removes_self_hits_and_full_length_aliases():
    results = {
        "snapgene:parent": _annotations(
            [
                {"db": "snapgene", "sseqid": "parent"},
                {"db": "snapgene", "sseqid": "alias"},
                {
                    "db": "snapgene",
                    "sseqid": "child",
                    "name": "T7 promoter",
                    "type": "promoter",
                    "qstart": 25,
                    "qend": 44,
                    "length": 19,
                    "slen": 19,
                },
            ]
        )
    }
    parents = {
        "snapgene:parent": {
            "db": "snapgene",
            "sseqid": "parent",
            "name": "parent",
            "type": "intron",
            "length": 100,
        }
    }

    report, missing = nested_rows(
        results,
        parents,
        {"snapgene": "blastn"},
    )

    assert report["nested_sseqid"].tolist() == ["child"]
    assert missing == []


def test_nested_rows_discards_nucleotide_hits_in_back_translated_fpbase():
    results = {
        "fpbase:parent": _annotations(
            [
                {"db": "fpbase", "sseqid": "parent"},
                {
                    "db": "snapgene",
                    "sseqid": "synthetic_dna_hit",
                    "qstart": 10,
                    "qend": 30,
                    "length": 20,
                },
                {
                    "db": "swissprot",
                    "sseqid": "protein_hit",
                    "qstart": 40,
                    "qend": 70,
                    "length": 30,
                },
            ]
        )
    }
    parents = {
        "fpbase:parent": {
            "db": "fpbase",
            "sseqid": "parent",
            "name": "parent",
            "type": "CDS",
            "length": 100,
        }
    }

    report, missing = nested_rows(
        results,
        parents,
        {"fpbase": "diamond", "snapgene": "blastn", "swissprot": "diamond"},
    )

    assert report["nested_sseqid"].tolist() == ["protein_hit"]
    assert missing == []


def test_figure_records_are_linear_and_convert_coordinates():
    report = pd.DataFrame(
        [
            {
                "parent_db": "snapgene",
                "parent_sseqid": "CMV_intron_(3)",
                "parent_name": "CMV intron",
                "parent_type": "intron",
                "parent_length": 110,
                "nested_db": "snapgene",
                "nested_sseqid": "T7_promoter",
                "nested_name": "T7 promoter",
                "nested_type": "promoter",
                "nested_strand": -1,
                "nested_start": 55,
                "nested_end": 74,
                "percent_identity": 100.0,
                "percent_match": 100.0,
                "evalue": 1e-20,
                "fragment": False,
            }
        ]
    )

    [record] = build_records(report)

    assert record["parsed"]["topology"] == "linear"
    assert record["parsed"]["length"] == 110
    assert record["children"][0]["status"] == "good"
    assert record["children"][0]["action"] == "keep"
    assert record["children"][0]["evalue"] == 1e-20
    assert record["parsed"]["features"][1] == {
        "type": "promoter",
        "start": 56,
        "end": 74,
        "strand": -1,
        "spansOrigin": False,
        "qualifiers": {
            "label": "T7 promoter",
            "database": "snapgene",
            "identity": "100.0",
            "match_length": "100.0",
            "fragment": "False",
            "note": "snapgene:T7_promoter",
        },
    }


def test_curated_pair_wins_before_fragment_suppression():
    decision = classify_row(
        {
            "parent_db": "snapgene",
            "parent_sseqid": "T5_promoter",
            "nested_db": "snapgene",
            "nested_sseqid": "lac_operator_(2)",
            "fragment": True,
        }
    )

    assert decision.status == "good"
    assert decision.action == "replace_parent"


def test_packaged_pair_overrides_expand_to_unique_source_pinned_keys():
    decisions = curated_decisions()

    assert all(all(key) for key in decisions)
    origin_parents = {
        parent
        for parent_db, parent, child_db, child in decisions
        if (parent_db, child_db, child) == ("snapgene", "Rfam", "RF00106")
    }
    assert origin_parents == {"ori", "p15A_ori", "ColA_ori", "CloDF13_ori", "RSF_ori"}
    assert (
        decisions[("snapgene", "OpIE_2_promoter", "swissprot", "P41708")].action
        == "suppress_child"
    )


def test_fragment_suppression_and_whole_child_keep_rules():
    base = {
        "parent_db": "snapgene",
        "parent_sseqid": "parent",
        "parent_type": "promoter",
        "nested_db": "snapgene",
        "nested_sseqid": "child",
        "nested_type": "terminator",
    }

    short_fragment = base | {
        "fragment": True,
        "nested_length": 60,
        "percent_identity": 100.0,
    }
    assert classify_row(short_fragment).action == "suppress_child"
    whole = classify_row(base | {"fragment": False})
    assert whole.status == "good"
    assert whole.action == "keep"
    assert whole.source == "rule:whole_child_keep"
    assert (
        classify_row(
            base
            | {
                "nested_type": "promoter",
                "fragment": True,
                "nested_length": 120,
                "percent_identity": 100.0,
            }
        ).status
        == "review"
    )


def test_high_confidence_cds_fragment_survives_as_provenance():
    decision = classify_row(
        {
            "parent_db": "snapgene",
            "parent_sseqid": "viral_promoter",
            "parent_type": "promoter",
            "nested_db": "swissprot",
            "nested_sseqid": "viral_protein",
            "nested_type": "CDS",
            "nested_length": 144,
            "percent_identity": 99.0,
            "fragment": True,
        }
    )

    assert decision.status == "good"
    assert decision.action == "keep"
    assert decision.source == "rule:high_confidence_cds_fragment"


def test_weak_or_short_cds_fragments_do_not_pass_general_keep_rule():
    base = {
        "parent_db": "snapgene",
        "parent_sseqid": "viral_promoter",
        "parent_type": "promoter",
        "nested_db": "swissprot",
        "nested_sseqid": "viral_protein",
        "nested_type": "CDS",
        "fragment": True,
    }

    assert (
        classify_row(base | {"nested_length": 87, "percent_identity": 100.0}).status
        == "bad"
    )
    assert (
        classify_row(base | {"nested_length": 180, "percent_identity": 94.9}).status
        == "bad"
    )


def test_reviewed_low_identity_viral_pairs_override_general_threshold():
    mmlv = classify_row(
        {
            "parent_db": "snapgene",
            "parent_sseqid": "MMLV_Psi",
            "parent_type": "misc_feature",
            "nested_db": "swissprot",
            "nested_sseqid": "P0DOG8",
            "nested_type": "CDS",
            "nested_length": 207,
            "percent_identity": 72.5,
            "fragment": True,
        }
    )
    opie2 = classify_row(
        {
            "parent_db": "snapgene",
            "parent_sseqid": "OpIE_2_promoter",
            "parent_type": "promoter",
            "nested_db": "swissprot",
            "nested_sseqid": "P41708",
            "nested_type": "CDS",
            "nested_length": 168,
            "percent_identity": 58.9,
            "fragment": True,
        }
    )

    assert mmlv.status == "good"
    assert mmlv.action == "replace_child"
    assert opie2.status == "bad"
    assert opie2.source == "curated_pair:opie2_ac152_wrong_frame"


def test_rfam_ncrna_and_compound_gene_components_are_kept():
    ncrna = classify_row(
        {
            "parent_db": "snapgene",
            "parent_sseqid": "viral_ltr",
            "parent_type": "LTR",
            "nested_db": "custom_cm",
            "nested_sseqid": "RF00001",
            "nested_type": "ncRNA",
            "fragment": True,
        }
    )
    component = classify_row(
        {
            "parent_db": "snapgene",
            "parent_sseqid": "compound_marker",
            "parent_type": "gene",
            "nested_db": "snapgene",
            "nested_sseqid": "marker_promoter",
            "nested_type": "promoter",
            "fragment": False,
        }
    )

    assert ncrna.source == "rule:ncrna_keep"
    assert component.source == "rule:compound_gene_component"
    assert ncrna.status == component.status == "good"


def test_boundary_extension_requires_edge_overhang_length_identity_and_evalue():
    base = {
        "parent_db": "snapgene",
        "parent_sseqid": "short_parent",
        "parent_type": "CDS",
        "parent_length": 300,
        "nested_db": "snapgene",
        "nested_sseqid": "long_child",
        "nested_type": "gene",
        "nested_start": 1,
        "nested_end": 298,
        "nested_length": 297,
        "nested_subject_length": 900,
        "nested_subject_start": 1,
        "nested_subject_end": 297,
        "percent_identity": 99.0,
        "evalue": 1e-40,
        "fragment": True,
    }

    assert classify_row(base).source == "rule:boundary_extension_keep"
    assert classify_row(base | {"nested_start": 20, "nested_end": 280}).status == "bad"
    assert classify_row(base | {"nested_subject_length": 320}).status == "bad"
    assert classify_row(base | {"evalue": 1e-3}).status == "bad"


def test_boundary_extension_respects_reverse_strand_subject_tail():
    base = {
        "parent_db": "snapgene",
        "parent_sseqid": "parent",
        "parent_type": "promoter",
        "parent_length": 300,
        "nested_db": "swissprot",
        "nested_sseqid": "child",
        "nested_type": "CDS",
        "nested_strand": -1,
        "nested_start": 0,
        "nested_end": 180,
        "nested_length": 180,
        "nested_subject_length": 600,
        "nested_subject_start": 1,
        "nested_subject_end": 180,
        "percent_identity": 99.0,
        "percent_match": 30.0,
        "evalue": 1e-40,
        "fragment": True,
    }

    # On the reverse query strand, the subject suffix continues through the
    # parent's left edge. The same subject span on the plus strand does not.
    assert classify_row(base).source == "rule:boundary_extension_keep"
    assert classify_row(base | {"nested_strand": 1}).source != (
        "rule:boundary_extension_keep"
    )


def test_near_complete_fragment_uses_coverage_with_evalue_as_support():
    base = {
        "parent_db": "snapgene",
        "parent_sseqid": "parent",
        "parent_type": "promoter",
        "nested_db": "swissprot",
        "nested_sseqid": "child",
        "nested_type": "CDS",
        "nested_length": 240,
        "percent_match": 80.0,
        "percent_identity": 84.0,
        "evalue": 1e-30,
        "fragment": True,
    }

    assert classify_row(base).source == "rule:near_complete_fragment_keep"
    assert classify_row(base | {"evalue": 1e-3}).status == "bad"
    # Exact same-source DNA does not need an arbitrary e-value cutoff to establish
    # that almost all of the child record is physically present.
    same_source = base | {
        "nested_db": "snapgene",
        "percent_identity": 100.0,
        "evalue": 1e-3,
    }
    assert classify_row(same_source).source == "rule:near_complete_fragment_keep"


def test_exact_short_non_cds_elements_survive_at_moderate_coverage():
    decision = classify_row(
        {
            "parent_db": "snapgene",
            "parent_sseqid": "2u_ori_(2)",
            "parent_type": "rep_origin",
            "parent_length": 983,
            "nested_db": "snapgene",
            "nested_sseqid": "FRT_(minimal)",
            "nested_type": "protein_bind",
            "nested_length": 20,
            "nested_subject_length": 34,
            "nested_subject_start": 15,
            "nested_subject_end": 34,
            "percent_identity": 100.0,
            "percent_match": 58.8,
            "evalue": 7.36e-4,
            "fragment": True,
        }
    )

    assert decision.status == "good"
    assert decision.source == "rule:exact_short_element_keep"
    below_threshold = {
        "parent_db": "snapgene",
        "parent_sseqid": "parent",
        "parent_type": "promoter",
        "nested_db": "snapgene",
        "nested_sseqid": "long_fragment",
        "nested_type": "terminator",
        "percent_identity": 100.0,
        "percent_match": 29.9,
        "fragment": True,
    }
    assert classify_row(below_threshold).action == "suppress_child"


def test_composite_component_only_fragment_is_suppressed_before_general_keep():
    row = {
        "parent_db": "snapgene",
        "parent_sseqid": "lac_operator",
        "parent_type": "protein_bind",
        "parent_length": 25,
        "nested_db": "snapgene",
        "nested_sseqid": "T5_promoter",
        "nested_type": "promoter",
        "nested_length": 19,
        "nested_subject_length": 45,
        "nested_subject_start": 21,
        "nested_subject_end": 39,
        "percent_identity": 100.0,
        "percent_match": 42.2,
        "evalue": 3.92e-5,
        "fragment": True,
    }

    decision = classify_row(row)

    assert decision.status == "bad"
    assert decision.action == "suppress_child"
    assert decision.source == "curated_region:t5_laco_component"


def test_composite_region_filter_is_unary_direction_agnostic_and_fail_open():
    t5_fragment = {
        "db": "snapgene",
        "sseqid": "T5_promoter",
        "sstart": 21,
        "send": 39,
        "fragment": True,
    }
    unaffected = [
        # A full record is never removed.
        t5_fragment | {"sstart": 1, "send": 45, "fragment": False},
        # Twelve distinctive T5 bases remain outside the embedded lacO interval.
        t5_fragment | {"sstart": 9, "send": 22},
        # A short interval wholly outside every component is informative, not covered.
        t5_fragment | {"sstart": 38, "send": 40},
        # Curated ids are exact and source-pinned.
        t5_fragment | {"db": "custom"},
    ]

    filtered = suppress_uninformative_composite_fragments(
        pd.DataFrame(
            [t5_fragment, t5_fragment | {"sstart": 39, "send": 21}, *unaffected]
        )
    )

    assert filtered.to_dict("records") == unaffected


def test_curated_fragment_region_suppresses_only_the_weak_pena_artifact():
    pena_fragment = {
        "db": "swissprot",
        "sseqid": "Q02940",
        "slen": 939,
        "sstart": 580,
        "send": 675,
        "pident": 68.8,
        "fragment": True,
    }
    unaffected = [
        pena_fragment | {"fragment": False},
        pena_fragment | {"pident": 70.1},
        pena_fragment | {"slen": 936},
        pena_fragment | {"sstart": 400, "send": 495},
        # Only 89.1% of this shifted interval overlaps the curated region.
        pena_fragment | {"sstart": 569, "send": 669},
        # A much larger hit containing the region is not the curated artifact.
        pena_fragment | {"sstart": 100, "send": 900},
        pena_fragment | {"db": "custom"},
    ]

    filtered = suppress_curated_fragment_artifacts(
        pd.DataFrame(
            [
                pena_fragment,
                pena_fragment | {"sstart": 675, "send": 580},
                # Exactly 90% of this interval overlaps the curated region.
                pena_fragment | {"sstart": 570, "send": 669},
                *unaffected,
            ]
        )
    )

    assert filtered.to_dict("records") == unaffected


def test_curated_fragment_region_is_visible_in_report_classification():
    decision = classify_row(
        {
            "parent_db": "snapgene",
            "parent_sseqid": "lac_promoter",
            "nested_db": "swissprot",
            "nested_sseqid": "Q02940",
            "nested_subject_length": 939,
            "nested_subject_start": 580,
            "nested_subject_end": 675,
            "percent_identity": 68.8,
            "fragment": True,
        }
    )

    assert decision.status == "bad"
    assert decision.action == "suppress_child"
    assert decision.source == "curated_fragment:puc_lac_region_penA"


def test_production_nested_suppression_is_conservative_and_fail_open():
    parent = {
        "db": "snapgene",
        "sseqid": "container",
        "type": "promoter",
        "qstart": 100,
        "qend": 500,
        "qlen": 1000,
        "score": 400.0,
        "fragment": False,
        "sframe": 1,
        "length": 400,
        "slen": 400,
        "sstart": 1,
        "send": 400,
        "pident": 100.0,
        "abs percmatch": 100.0,
        "evalue": 0.0,
    }
    weak_child = parent | {
        "db": "snapgene",
        "sseqid": "short_child",
        "type": "terminator",
        "qstart": 200,
        "qend": 240,
        "score": 40.0,
        "fragment": True,
        "length": 40,
        "slen": 200,
        "sstart": 1,
        "send": 40,
        "abs percmatch": 20.0,
    }
    near_complete = weak_child | {
        "sseqid": "near_complete",
        "qstart": 300,
        "qend": 380,
        "length": 80,
        "slen": 100,
        "send": 80,
        "abs percmatch": 80.0,
    }

    filtered = suppress_nested_fragments(
        pd.DataFrame([parent, weak_child, near_complete])
    )

    assert filtered["sseqid"].tolist() == ["container", "near_complete"]

    # A second containing parent of the same kind as the child yields review; any
    # non-suppression decision wins rather than allowing a broad parent to erase it.
    reviewing_parent = parent | {
        "sseqid": "same_kind_container",
        "type": "terminator",
        "qstart": 150,
        "qend": 260,
    }
    kept = suppress_nested_fragments(
        pd.DataFrame([parent, reviewing_parent, weak_child])
    )
    assert "short_child" in set(kept["sseqid"])


def test_production_suppression_handles_origin_spanning_containment():
    parent = {
        "db": "snapgene",
        "sseqid": "wrapped_parent",
        "type": "promoter",
        "qstart": 900,
        "qend": 100,
        "qlen": 1000,
        "score": 500.0,
        "fragment": False,
        "sframe": 1,
        "length": 200,
        "slen": 200,
        "sstart": 1,
        "send": 200,
        "pident": 100.0,
        "abs percmatch": 100.0,
        "evalue": 0.0,
    }
    child = parent | {
        "sseqid": "weak_child",
        "type": "terminator",
        "qstart": 950,
        "qend": 980,
        "score": 30.0,
        "fragment": True,
        "length": 30,
        "slen": 300,
        "sstart": 101,
        "send": 130,
        "abs percmatch": 10.0,
    }

    assert suppress_nested_fragments(pd.DataFrame([parent, child]))[
        "sseqid"
    ].tolist() == ["wrapped_parent"]
    across_origin = child | {"qstart": 970, "qend": 30}
    assert suppress_nested_fragments(pd.DataFrame([parent, across_origin]))[
        "sseqid"
    ].tolist() == ["wrapped_parent"]


def test_curated_suppress_child_applies_in_production_geometry():
    parent = {
        "db": "snapgene",
        "sseqid": "OpIE_2_promoter",
        "type": "promoter",
        "qstart": 100,
        "qend": 400,
        "qlen": 1000,
        "score": 400.0,
        "fragment": False,
        "sframe": 1,
        "length": 300,
        "slen": 300,
        "sstart": 1,
        "send": 300,
        "pident": 100.0,
        "abs percmatch": 100.0,
        "evalue": 0.0,
    }
    child = parent | {
        "db": "swissprot",
        "sseqid": "P41708",
        "type": "CDS",
        "qstart": 160,
        "qend": 328,
        "score": 100.0,
        "fragment": True,
        "length": 168,
        "slen": 900,
        "sstart": 1,
        "send": 168,
        "pident": 58.9,
        "abs percmatch": 18.7,
        "evalue": 1e-20,
    }

    assert suppress_nested_fragments(pd.DataFrame([parent, child]))[
        "sseqid"
    ].tolist() == ["OpIE_2_promoter"]


def test_broken_fpbase_parents_are_source_corrections():
    for parent, child in (("stayrose", "staygold"), ("wiphy2", "wi-phy")):
        decision = classify_row(
            {
                "parent_db": "fpbase",
                "parent_sseqid": parent,
                "nested_db": "fpbase",
                "nested_sseqid": child,
                "fragment": False,
            }
        )
        assert decision.status == "bad"
        assert decision.action == "replace_parent"


def test_gcamp_calmodulin_fragment_is_kept_as_fusion_provenance():
    decision = classify_row(
        {
            "parent_db": "fpbase",
            "parent_sseqid": "gcamp2",
            "nested_db": "swissprot",
            "nested_sseqid": "P62146",
            "fragment": True,
        }
    )

    assert decision.status == "good"
    assert decision.action == "relabel_child"
    assert decision.source == "curated_pair:gcamp2_calmodulin_component"


def test_current_audit_is_partitioned_by_policy():
    report = pd.read_csv(ROOT / "docs" / "nested-feature-audit.csv")
    classified = classify_report(report)

    assert len(classified) == 478
    assert classified["rule_status"].value_counts().to_dict() == {
        "bad": 263,
        "good": 215,
    }
