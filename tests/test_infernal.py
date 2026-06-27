"""Tests for Infernal output parsing and coordinate handling."""

import pandas as pd

from plannotate._filter import _normalize_coordinates
from plannotate._tools.infernal import parse_output

INFERNAL_COLUMNS = [
    "#idx",
    "target name",
    "accession",
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
