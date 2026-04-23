# pLannotate Repo Review

## Top Changes

1. Decouple the annotation engine from Streamlit.

   `plannotate/annotate.py` imports `streamlit` in core logic, uses `st.cache`, `st.progress`, and `st.error` inside library paths. Make `annotate()` pure and let the Streamlit app pass an optional progress callback/cache wrapper.

2. Harden subprocess execution.

   `plannotate/annotate.py` runs BLAST/DIAMOND/Infernal without `check=True`, timeout handling, or stderr surfacing. It also uses `shell=True` with `rg`. Centralize command execution with explicit argv lists, timeouts, dependency checks, and clear user-facing errors.

3. Replace the base-by-base overlap filter.

   `plannotate/annotate.py` builds a dataframe with one column per base position, then scans it for overlap resolution. That is fragile and can scale badly, especially because CLI batch bypasses max length. Rewrite this around interval operations.

4. Make tests fast, separated, and observable.

   The suite mixes unit tests with external database/tool integration tests. A local `python -m pytest tests/test_units.py -q` run hung silently; even `SIGKILL` left PID `2568` in uninterruptible wait. Add integration markers, per-test timeouts, fake subprocess fixtures for unit tests, and keep external DB tests in a separate CI job.

   Status: addressed by adding an integration marker gate, per-test timeout guard, lazy imports for deselected integration tests, and separate unit/integration CI steps.

5. Make database installation/versioning robust.

   `resources.databases_exist()` only checks whether `data/BLAST_dbs/` exists. `download_databases()` downloads into the installed package directory and removes files via shell tools. Use a user cache dir, checksum/manifest validation, atomic extraction, and verify required DB files.

6. Modernize the Streamlit/Bokeh dependency path.

   Dependencies are tightly pinned around older UI APIs, and the CLI calls private Streamlit API `_main_run`. Upgrade deliberately and remove private/deprecated calls.

7. Tighten CLI/file handling.

   The CLI builds paths with string interpolation and uses `print`/`sys.exit` for errors. Use `pathlib`, create/validate output dirs, use Click exceptions, and fix help text typos.

8. Clean up plotting internals after core stability.

   `plannotate/bokeh_plot.py` has a "weird error" fallback for levels, mutates dataframes, and uses older Bokeh parameter names. This is lower priority than annotation/test robustness, but worth doing.

## Verification Notes

- `python -m ruff check .` passed.
- Full pytest was not completed because of the stuck pytest process noted above.
