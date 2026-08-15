"""Tests for Infernal output parsing and coordinate handling."""

import pandas as pd

from plannotate._filter import _normalize_coordinates
from plannotate._tools.infernal import model_lengths, parse_alignments, parse_output

# One cmscan alignment block per hit, trimmed to the lines the parser reads. The
# second block is wrapped, which is how any hit longer than the report width arrives.
ALIGNMENT_REPORT = """\
Query:       q0  [L=12030]
>> SAM_riboswitch  SAM riboswitch
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (1) !   7.6e-19   90.4   0.0  cm        3       94 []           1          92 + .. 0.97    no 0.65

                                                                    NC
           :::::<<<<<<<____>>>>>>>::: CS
  SAM       3 aguAUUUGGuggCuGcGcuCuucua 27
             AGUAUUUGGU::CUGCGCUCU CU+
     q0    1 AGUAUUUGGUAUCUGCGCUCUGCUG 25
             ************************* PP

>> Reverse_hit  reverse strand hit
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (2) !   2.4e-05   31.0   0.0  cm       30       10 []          90          70 - .. 0.88    no 0.48

           <<<<<____ CS
  Rev      30 aguAUUUGG 22
             AGUAUUUGG
     q0   90 AGUAUUUGG 82
             ********* PP

           >>>>> CS
  Rev      21 uggCu 10
             UGGCU
     q0   81 UGGCU 70
             ***** PP
"""

INFERNAL_COLUMNS = [
    "#idx",
    "target name",
    "accession",
    "query name",
    "clan name",
    "seq from",
    "seq to",
    "mdl from",
    "mdl to",
    "strand",
    "score",
    "E-value",
    "description of target",
]


def _write_infernal_fixture(tmp_path, rows):
    widths = []
    for index, column in enumerate(INFERNAL_COLUMNS):
        row_width = max((len(str(row[index])) for row in rows), default=0)
        widths.append(max(len(column), row_width, 8) + 2)

    header = "".join(
        column.ljust(width)
        for column, width in zip(INFERNAL_COLUMNS, widths, strict=True)
    )
    dividers = []
    for index, width in enumerate(widths):
        marker = "#" + "-" * (width - 2) if index == 0 else "-" * (width - 1)
        dividers.append(marker.ljust(width))

    lines = [header, "".join(dividers)]
    lines.extend(
        "".join(
            str(value).ljust(width) for value, width in zip(row, widths, strict=True)
        )
        for row in rows
    )
    path = tmp_path / "infernal.tbl"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_parse_infernal_empty_tblout(tmp_path):
    parsed = parse_output(_write_infernal_fixture(tmp_path, []))

    assert parsed.empty
    assert "accession" not in parsed.columns
    assert "clan name" not in parsed.columns


def test_parse_infernal_preserves_one_based_inclusive_coordinates(tmp_path):
    tblout = _write_infernal_fixture(
        tmp_path,
        [
            [
                1,
                "SAM_riboswitch",
                "RF00162",
                "query",
                "CL00001",
                1,
                92,
                3,
                94,
                "+",
                "42.5",
                "1e-12",
                "SAM riboswitch",
            ],
            [
                2,
                "Reverse_hit",
                "-",
                "query",
                "-",
                90,
                70,
                30,
                10,
                "-",
                "31.0",
                "2e-05",
                "reverse strand hit",
            ],
        ],
    )

    parsed = parse_output(tblout)

    assert parsed.loc[0, "name"] == "SAM riboswitch"
    assert parsed.loc[0, "qstart"] == 1
    assert parsed.loc[0, "qend"] == 92
    assert parsed.loc[0, "length"] == 92
    assert parsed.loc[0, "sframe"] == 1
    assert parsed.loc[1, "qstart"] == 70
    assert parsed.loc[1, "qend"] == 90
    assert parsed.loc[1, "sframe"] == -1
    # sseqid is the stable Rfam accession, with a fallback to the model name when
    # the model carries no accession (the "-" field in the second fixture row)
    assert parsed.loc[0, "sseqid"] == "RF00162"
    assert parsed.loc[1, "sseqid"] == "Reverse hit"
    assert "#idx" not in parsed.columns


def test_shared_coordinate_normalizer_converts_infernal_once():
    hits = pd.DataFrame({"qstart": [1], "qend": [92]})

    normalized = _normalize_coordinates(hits)

    assert normalized.loc[0, "qstart"] == 0
    assert normalized.loc[0, "qend"] == 91


def test_parse_alignments_reads_structure_and_accuracy(tmp_path):
    path = tmp_path / "alignments.txt"
    path.write_text(ALIGNMENT_REPORT)

    alignments = parse_alignments(path)

    # keyed on the query and the raw coordinates cmscan prints, which run descending on
    # the minus strand -- the table parser sorts them only afterwards
    assert set(alignments) == {
        ("q0", "SAM_riboswitch", "1", "92"),
        ("q0", "Reverse_hit", "90", "70"),
    }
    structure, accuracy = alignments[("q0", "SAM_riboswitch", "1", "92")]
    assert structure == ":::::<<<<<<<____>>>>>>>:::"
    assert accuracy == 0.97
    # a wrapped hit is stitched back into one string rather than truncated at the
    # first block
    assert alignments[("q0", "Reverse_hit", "90", "70")][0] == "<<<<<____>>>>>"


def test_parse_alignments_degrades_instead_of_raising(tmp_path):
    # the alignment report is not a documented interface, so a shape the parser does
    # not recognise must cost the structure and nothing else
    garbled = tmp_path / "garbled.txt"
    garbled.write_text(">> model  description\n  (1) ! not a number\n")

    assert parse_alignments(garbled) == {}
    assert parse_alignments(tmp_path / "absent.txt") == {}


def test_parse_alignments_separates_queries_hitting_a_model_alike(tmp_path):
    # a batched search puts every sequence through one cmscan, and two plasmids sharing
    # a backbone can hit the same model at the same offset. Keying on the model and
    # coordinates alone would let the second block overwrite the first.
    path = tmp_path / "alignments.txt"
    path.write_text(
        "Query:       q0  [L=100]\n"
        ">> RNAI  RNAI\n"
        "  (1) !   1e-12   73.7   0.0  cm        1      4 []           1           4 + .. 0.93    no 0.48\n"
        "           <<>> CS\n"
        "Query:       q1  [L=100]\n"
        ">> RNAI  RNAI\n"
        "  (1) !   1e-12   73.7   0.0  cm        1      4 []           1           4 + .. 0.72    no 0.48\n"
        "           :::: CS\n"
    )

    alignments = parse_alignments(path)

    assert alignments[("q0", "RNAI", "1", "4")] == ("<<>>", 0.93)
    assert alignments[("q1", "RNAI", "1", "4")] == ("::::", 0.72)


def test_one_unreadable_block_does_not_discard_the_others(tmp_path):
    # parsing is per hit, so a single block in a shape the parser does not recognise
    # costs that hit's structure and leaves every other hit intact
    path = tmp_path / "alignments.txt"
    path.write_text(ALIGNMENT_REPORT + "\n>> Broken  truncated row\n  (1) ! oops\n")

    alignments = parse_alignments(path)

    assert len(alignments) == 2
    assert alignments[("q0", "SAM_riboswitch", "1", "92")][1] == 0.97


def test_parse_output_attaches_structure_identity_and_model_length(tmp_path):
    tblout = _write_infernal_fixture(
        tmp_path,
        [
            [
                1,
                "SAM_riboswitch",
                "RF00162",
                "q0",
                "CL00001",
                1,
                92,
                3,
                94,
                "+",
                "42.5",
                "1e-12",
                "SAM riboswitch",
            ],
            [
                2,
                "Reverse_hit",
                "-",
                "q0",
                "-",
                90,
                70,
                30,
                10,
                "-",
                "31.0",
                "2e-05",
                "reverse strand hit",
            ],
        ],
    )
    alignments = tmp_path / "alignments.txt"
    alignments.write_text(ALIGNMENT_REPORT)
    cm = tmp_path / "models.cm"
    cm.write_text("NAME     SAM_riboswitch\nACC      RF00162\nCLEN     108\n")

    parsed = parse_output(tblout, alignments, cm)

    assert parsed.loc[0, "structure"] == ":::::<<<<<<<____>>>>>>>:::"
    # identity carries the alignment's posterior-probability accuracy, because a
    # covariance model scores structure rather than base identity
    assert parsed.loc[0, "pident"] == 97.0
    # slen is the model's full width, so the partial hit (model 3-94 of 108) is not
    # reported as a full-length match
    assert parsed.loc[0, "slen"] == 108
    # a model absent from the CM file falls back to the aligned span
    assert parsed.loc[1, "slen"] == 21


def test_parse_output_without_an_alignment_report_keeps_working(tmp_path):
    tblout = _write_infernal_fixture(
        tmp_path,
        [
            [
                1,
                "SAM_riboswitch",
                "RF00162",
                "q0",
                "CL00001",
                1,
                92,
                3,
                94,
                "+",
                "42.5",
                "1e-12",
                "SAM riboswitch",
            ]
        ],
    )

    parsed = parse_output(tblout)

    assert parsed.loc[0, "structure"] == ""
    assert parsed.loc[0, "pident"] == 100.0


def test_model_lengths_ignores_the_embedded_hmm_filter(tmp_path):
    # a pressed database follows each covariance model with the HMM filter built from
    # it, repeating NAME and ACC but reporting LENG instead of CLEN. Carrying those
    # repeats forward would bind every model to the *next* model's length.
    cm = tmp_path / "models.cm"
    cm.write_text(
        "NAME     RNAI\nACC      RF00106\nCLEN     102\n"
        "HMMER3/f\nNAME  RNAI\nACC   RF00106\nLENG  102\n"
        "NAME     5S_rRNA\nACC      RF00001\nCLEN     120\n"
        "HMMER3/f\nNAME  5S_rRNA\nACC   RF00001\nLENG  120\n"
    )

    lengths = model_lengths(cm)

    assert lengths["RF00106"] == 102
    assert lengths["RNAI"] == 102
    assert lengths["RF00001"] == 120
    assert lengths["5S_rRNA"] == 120


def test_model_lengths_result_cannot_be_corrupted_by_a_caller(tmp_path):
    cm = tmp_path / "models.cm"
    cm.write_text("NAME     RNAI\nACC      RF00106\nCLEN     102\n")

    model_lengths(cm).clear()

    assert model_lengths(cm)["RF00106"] == 102
