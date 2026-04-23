from plannotate.infernal import parse_infernal


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
        row_width = max(len(str(row[index])) for row in rows) if rows else 0
        widths.append(max(len(column), row_width, 8) + 2)

    header = "".join(
        column.ljust(width) for column, width in zip(INFERNAL_COLUMNS, widths)
    )
    dividers = []
    for index, width in enumerate(widths):
        if index == 0:
            dividers.append(("#" + "-" * (width - 2)).ljust(width))
        else:
            dividers.append(("-" * (width - 1)).ljust(width))
    divider = "".join(dividers)

    lines = [header, divider]
    for row in rows:
        lines.append(
            "".join(str(value).ljust(width) for value, width in zip(row, widths))
        )

    path = tmp_path / "infernal.tbl"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_parse_infernal_empty_tblout(tmp_path):
    tblout = _write_infernal_fixture(tmp_path, [])

    parsed = parse_infernal(tblout)

    assert parsed.empty
    assert "accession" not in parsed.columns
    assert "clan name" not in parsed.columns


def test_parse_infernal_normalizes_rows(tmp_path):
    tblout = _write_infernal_fixture(
        tmp_path,
        [
            [
                1,
                "SAM_riboswitch",
                "RF00162",
                "CL00001",
                5,
                25,
                3,
                23,
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

    parsed = parse_infernal(tblout)

    assert parsed.loc[0, "sseqid"] == 1
    assert parsed.loc[0, "Feature"] == "SAM riboswitch"
    assert parsed.loc[0, "Description"] == "Accession: RF00162 - SAM riboswitch"
    assert parsed.loc[0, "qstart"] == 5
    assert parsed.loc[0, "qend"] == 25
    assert parsed.loc[0, "sframe"] == 1
    assert parsed.loc[0, "length"] == 21
    assert parsed.loc[0, "slen"] == 21
    assert parsed.loc[0, "pident"] == 100
    assert parsed.loc[0, "qseq"] == ""

    assert parsed.loc[1, "Feature"] == "Reverse hit"
    assert parsed.loc[1, "Description"] == "Accession:   - reverse strand hit"
    assert parsed.loc[1, "qstart"] == 70
    assert parsed.loc[1, "qend"] == 90
    assert parsed.loc[1, "sframe"] == -1
    assert parsed.loc[1, "length"] == 21
    assert parsed.loc[1, "slen"] == 21
