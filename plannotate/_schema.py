"""Internal column definitions for annotation data frames."""

ADAPTER_COLUMNS = [
    "sseqid",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sframe",
    "evalue",
    "qseq",
    "length",
    "slen",
    "pident",
    "qlen",
]

BLAST_COLUMNS = [
    "sseqid",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sframe",
    "score",
    "evalue",
    "qseq",
    "length",
    "slen",
    "pident",
    "qlen",
    "db",
]

FEATURE_COLUMNS = [
    "name",
    "blurb",
    "type",
    "priority",
    "percmatch",
    "abs percmatch",
    "pi_permatch",
    "fragment",
]

PROCESSING_COLUMNS = [
    "wiggle",
    "wstart",
    "wend",
    "kind",
    "qstart_dup",
    "qend_dup",
]

ANNOTATION_COLUMNS = BLAST_COLUMNS + FEATURE_COLUMNS + PROCESSING_COLUMNS

CSV_COLUMNS = [
    "sseqid",
    "qstart",
    "qend",
    "sframe",
    "pident",
    "slen",
    "length",
    "abs percmatch",
    "fragment",
    "db",
    "name",
    "type",
    "blurb",
    "qseq",
]

CSV_COLUMN_NAMES = {
    "qstart": "start location",
    "qend": "end location",
    "sframe": "strand",
    "pident": "percent identity",
    "slen": "full length of feature in db",
    "qseq": "sequence",
    "length": "length of found feature",
    "abs percmatch": "percent match length",
    "db": "database",
    "name": "Feature",
    "type": "Type",
    "blurb": "Description",
}
