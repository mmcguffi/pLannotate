# Annotation controls

See `DIFFERENCES.md` for the current refactor-versus-1.2.5 comparison organized
for manual inspection.

These CSV and GenBank files record the annotation behavior of the published
`plannotate==1.2.5` implementation (tag `v1.2.5`, commit `61ed152`) for every
FASTA distributed in `plannotate/data/fastas`.

The controls cover:

- regular mode for all 10 example plasmids;
- detailed mode for all 10 example plasmids;
- linear mode for all 10 example plasmids;
- detailed and linear mode together for `pXampl3`.

The controls were generated from a clean environment created directly from the
Bioconda `plannotate=1.2.5` package (`pyhdfd78af_0`), then by running its
unmodified `plannotate batch` command. The databases bundled in that package are
byte-identical to the canonical archive downloaded by 1.2.5's `setupdb` command
from the `v1.2.0` GitHub release and match the release SHA-256. Exact package,
database, and search-tool versions are recorded in `regression-context.json`.
`current-database-manifest.json` records the database bundle expected when
running the refactored implementation.

GenBank comparisons cover sequence, topology, feature locations, types, and
qualifiers; generated dates and version comments are intentionally ignored. CSV
comparisons cover every column and value except that numeric Rfam IDs such as
`5.0` and `5` are treated as equivalent.

Changes emit pytest warnings and are reported as expected failures (`XFAIL`) by
default. This makes them visible in test runners such as VS Code without blocking
normal development. Each changed test also prints the added and removed CSV rows
and GenBank features for inspection. To make control changes fail the test run,
set:

```bash
PLANNOTATE_STRICT_ANNOTATION_CONTROLS=1 pytest --run-integration
```

The equivalent command-line option is `--strict-annotation-controls`. The local
VS Code test configuration enables this option because VS Code does not reliably
display captured output for XFAIL results. CI and ordinary pytest runs leave it
disabled.

Vanilla 1.2.5 itself produces two duplicate, untyped FPbase `11.0` annotations
over the fluorescent-protein feature in several detailed-mode controls. This was
verified directly against its canonical database archive. The fixtures preserve
that historical output even where the refactored SQLite metadata and overlap
filtering now produce cleaner results.
