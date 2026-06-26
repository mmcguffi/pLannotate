# Architecture

`annotate` is the single annotation pipeline. It collects candidate hits from all
configured sources, enriches them with feature metadata, resolves conflicts, and
returns the final annotation table. External tools remain isolated behind small
adapters because their commands and output formats change independently.

## Annotation pipeline

1. `validation` reads and validates one FASTA or GenBank record.
2. `annotate` loads the configured annotation sources.
3. `_concurrency` divides the core budget in YAML order and runs sources concurrently.
4. `_tools/methods` selects an integration from its function mappings.
5. `_tools/blast`, `_tools/diamond`, or `_tools/infernal` runs one executable and
   parses its output.
6. `annotate` attaches descriptions and priority, then asks `_filter` to score
   candidates and resolve overlaps.
7. `annotate` finalizes feature classification, coordinates, and query orientation.
8. `models` converts rows to `Feature` and `Construct` objects.

`_tools` contains the integrations and their shared subprocess, temporary-file, and
tabular-output helpers.
`_concurrency` owns both core allocation and enforcement of the allocated thread count
in tool parameters.

## Adding an annotation source

### Another database using an existing tool

Add one entry to `data/data/databases.yml`. No Python dispatch or concurrency code
needs to change. The entry selects a registered `method`, its location and parameters,
its filtering priority, and where feature descriptions come from.

### A new external tool or caller

1. Add one focused integration such as `_tools/prodigal.py` or `_tools/isescan.py`.
2. Implement `search(sequence, config, threads) -> pandas.DataFrame` in that module.
   The adapter owns command construction, output parsing, and tool-specific coordinate
   conversion.
3. Return the columns in `_schema.ADAPTER_COLUMNS`: `sseqid`, `qstart`, `qend`,
   `sstart`, `send`, `sframe`, `evalue`, `qseq`, `length`, `slen`, `pident`, and
   `qlen`. If `details.location` is `None`, also return `name`, `type`, and `blurb`;
   otherwise return `sseqid` values that match the configured descriptions database.
4. Register the function and its optional database path convention in
   `_tools.methods`.
5. Add parser unit tests and one integration test for the executable.

Database-backed detectors and database-free callers use the same adapter boundary. A
caller that needs no database maps its method to `None` in `DATABASE_DIRECTORIES`;
configuration loading then adds no database paths.

Adapter coordinates are 1-based and inclusive, matching external search output. A
caller without alignment statistics can use `evalue=0`, `pident=100`, `sstart=1`,
`send=length`, and `slen=length`; shared filtering converts coordinates and calculates
the final score.

## Dependency rules

- Tool adapters do not load descriptions, rank hits, or render output.
- Concurrency does not know which tools exist.
- Runtime concurrency uses the standard library; Snakemake is a database-build extra.
- Cached DataFrames are never returned directly to mutating callers.
- Optional Bokeh imports stay behind plotting APIs.
- Importing `plannotate` does not configure application logging.

## Public API

The supported top-level API is `plannotate.Construct`, `plannotate.Feature`, and
`plannotate.build_databases`. Modules beginning with `_` are implementation details.
