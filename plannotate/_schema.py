"""Internal column definitions for annotation data frames."""

ADAPTER_COLUMNS = [
    "qseqid",
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

# The canonical annotation row, in order: search statistics, feature metadata, then
# internal processing columns. This list defines the column order of every full
# annotation DataFrame and the CSV/GenBank round-trip; ``Feature.to_dict`` is driven
# by it. To add a column, add it here (and map it to a Feature field in models.py).
ANNOTATION_COLUMNS = [
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
    "name",
    "blurb",
    "type",
    "priority",
    "percmatch",
    "abs percmatch",
    "pi_permatch",
    "wiggle",
    "wstart",
    "wend",
    "kind",
    "qstart_dup",
    "qend_dup",
    "fragment",
]

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
