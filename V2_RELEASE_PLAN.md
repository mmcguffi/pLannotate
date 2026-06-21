# pLannotate 2.0.0 Snapshot Release Plan

## Decision

`big-refactor` will not be rebased or merged into `master` with Git's three-way
merge machinery. The branch is a replacement architecture with 84 commits since
its common ancestor with `master`, and a trial merge reports conflicts across
core code, packaging, workflows, and tests.

Instead, pLannotate 2.0.0 will be published as one snapshot commit:

- The commit's parent will be the latest `master` commit.
- The commit's file tree will exactly match the completed and tested
  `big-refactor` tree.
- The detailed development history will remain available on `big-refactor`.
- No history on `big-refactor` or `master` will be rewritten.

This intentionally treats 2.0.0 as a replacement of the 1.x architecture. It
avoids artificial merge-conflict resolution while retaining the repository's
existing history.

## Work to port before the snapshot

The following post-divergence `master` work must be manually adapted to the new
architecture. The original commits are references and should not be
cherry-picked wholesale.

Implementation status on `big-refactor`:

- Packaging, Python compatibility, Ruff, mypy, and CI have been ported.
- The Rfam coordinate fix and its unit/integration regressions have been ported.
- The deterministic 2.x database bundle and release-upload workflow have been
  added.
- Full external-tool integration remains an acceptance gate and requires the
  canonical database bundle configured through `PLANNOTATE_DATABASE_URL` and
  `PLANNOTATE_DATABASE_SHA256`.

### `98d6f14` — Python compatibility and packaging

- Replace `setup.py` with a `pyproject.toml` appropriate for the refactored
  package and its current CLI entry point.
- Consolidate runtime, optional plotting, testing, and lint dependencies.
- Support Python 3.10 through the currently supported Python versions.
- Validate pandas, NumPy, Bokeh, and Streamlit compatibility instead of merely
  copying the 1.x dependency constraints.
- Update installation documentation and CI installation commands.

### `ef683ce` — Rfam coordinate correction

- Remove duplicate Infernal coordinate normalization. Coordinates currently
  pass through both `plannotate/infernal.py` and
  `plannotate/filter_annotations.py`.
- Adjust sequence slicing to preserve the complete matched sequence.
- Port the issue-60 regression cases for left-edge, padded, linear, circular,
  and GenBank round-trip behavior.

### `fa77888` — tests, typing, logging, Ruff, and mypy

- Port test behavior, not old module structure.
- Separate fast unit tests from database/tool-dependent integration tests.
- Add deterministic test timeouts and an explicit integration-test switch.
- Configure and run Ruff and mypy against the refactored modules.
- Retain the refactor's logging design and port only missing behavior.
- Add CI jobs for unit tests, integration tests, Ruff, and mypy.

### `32801d0` — release database assets

- Do not reuse the workflow unchanged: it uploads the old v1.2.0 database
  bundle.
- Define the canonical 2.0.0 database artifact produced by the new gathering
  pipeline.
- Build or obtain that artifact reproducibly, verify its checksum, and attach
  it to the matching release.
- Verify that runtime database download and extraction use the same asset name,
  layout, and versioning policy.

The later Python-constraint commits (`f355fe4` and `61ed152`) reduce to a final
minimum constraint of Python 3.10 with no speculative upper bound. Actual
support is determined by the tested CI matrix.

## Pre-snapshot acceptance gates

Before creating the snapshot, the tip of `big-refactor` must:

1. Report package version `2.0.0` from the canonical package metadata.
2. Install successfully from a clean environment.
3. Pass the supported Python-version unit-test matrix.
4. Pass integration tests with external tools and the 2.0.0 databases.
5. Pass Ruff and mypy.
6. Build wheel and source-distribution artifacts successfully.
7. Exercise CLI, Python API, plotting, GenBank, CSV, and notebook paths.
8. Produce and consume the intended release database artifact.
9. Have a clean working tree and a pushed recovery branch.

## Creating the snapshot commit

Use a separate, clean worktree so replacing the index and working tree cannot
disturb the development checkout. Replace `<final-refactor-commit>` with the
reviewed commit SHA rather than relying on a moving branch name.

```bash
git fetch upstream

git worktree add ../pLannotate-v2 \
  -b v2-release upstream/master

cd ../pLannotate-v2

git read-tree --reset -u <final-refactor-commit>

git status
git diff --cached --stat
git diff --cached --check

git commit -m "Release pLannotate 2.0.0"
```

`git read-tree --reset -u` does not merge. It stages the exact tree from the
specified refactor commit on top of the latest `master` parent. Consequently,
the resulting commit is a normal Git commit with one parent and no merge
conflicts.

Do not run the `read-tree` command in a dirty or valuable worktree.

## Post-snapshot verification

The snapshot tree should equal the approved refactor tree exactly:

```bash
test "$(git rev-parse HEAD^{tree})" = \
  "$(git rev-parse <final-refactor-commit>^{tree})"
```

Then repeat all acceptance gates on `v2-release`. Testing only
`big-refactor` is insufficient because the release branch is the artifact that
will be reviewed and merged.

Review the complete replacement diff against its parent:

```bash
git diff --stat HEAD^
git diff --check HEAD^
git status --short
```

Push the snapshot branch and open the 2.0.0 pull request:

```bash
git push -u upstream v2-release
```

After review and CI pass, merge the commit without squashing it again. Tag the
resulting `master` commit as `v2.0.0` and publish the release assets through the
verified release workflow.

## Recovery and audit trail

- `big-refactor` remains the authoritative detailed development history.
- Record `<final-refactor-commit>` in the 2.0.0 pull-request description.
- Keep the refactor branch until after the 2.0.0 release has been validated in
  production.
- If the snapshot needs changes, make them first on `big-refactor`, rerun the
  gates, and recreate the unmerged snapshot commit. Once review begins, normal
  follow-up commits on `v2-release` are also acceptable and preserve review
  context.
- The snapshot commit can be reverted normally because it has a single parent.
