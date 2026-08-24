#!/usr/bin/env python3
"""Build a self-contained genes-d3 viewer for the nested-feature audit."""

import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import Any, cast

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from plannotate._nested import as_bool, classify_report  # noqa: E402

DEFAULT_INPUT = ROOT / "docs" / "nested-feature-audit.csv"
DEFAULT_DECISIONS = ROOT / "docs" / "nested-feature-decisions.csv"
DEFAULT_OUTPUT = ROOT / "figures" / "nested-feature-audit-linear"

GENES_D3_ASSETS = (
    "css/genes-d3.css",
    "js/parser.js",
    "js/data.js",
    "js/plot-core.js",
    "js/plot-linear.js",
    "js/plot-circular.js",
    "js/plot-stylized.js",
    "js/plot.js",
    "js/export.js",
)


def build_records(report: pd.DataFrame) -> list[dict[str, Any]]:
    """Convert audit pairs into the parsed feature objects genes-d3 consumes."""
    if "rule_status" not in report.columns:
        report = classify_report(report)
    records: list[dict[str, Any]] = []
    parent_columns = [
        "parent_db",
        "parent_sseqid",
        "parent_name",
        "parent_type",
        "parent_length",
    ]
    for parent, children in report.groupby(parent_columns, sort=False):
        parent_values = cast(tuple[object, object, object, object, object], parent)
        database, sequence_id, name, feature_type = map(str, parent_values[:4])
        length = int(str(parent_values[4]))
        features = [
            {
                "type": str(feature_type),
                "start": 1,
                "end": length,
                "strand": 1,
                "spansOrigin": False,
                "qualifiers": {
                    "label": str(name),
                    "database": str(database),
                    "identity": "100.0",
                    "match_length": "100.0",
                    "fragment": "False",
                    "note": "Source feature under audit",
                },
            }
        ]
        child_rows: list[dict[str, Any]] = []
        for _, child in children.iterrows():
            is_fragment = as_bool(child["fragment"])
            child_row: dict[str, Any] = {
                "db": str(child["nested_db"]),
                "sseqid": str(child["nested_sseqid"]),
                "name": str(child["nested_name"]),
                "type": str(child["nested_type"]),
                "strand": int(child["nested_strand"]),
                "start": int(child["nested_start"]),
                "end": int(child["nested_end"]),
                "identity": float(child["percent_identity"]),
                "match": float(child["percent_match"]),
                "evalue": float(child.get("evalue", float("nan"))),
                "fragment": is_fragment,
                "status": str(child["rule_status"]),
                "action": str(child["suggested_action"]),
                "rationale": str(child["decision_rationale"]),
                "decision_source": str(child["decision_source"]),
            }
            child_rows.append(child_row)
            features.append(
                {
                    "type": child_row["type"],
                    # Audit coordinates are zero-based, half-open; genes-d3's
                    # parsed GenBank representation is one-based, closed.
                    "start": child_row["start"] + 1,
                    "end": child_row["end"],
                    "strand": child_row["strand"],
                    "spansOrigin": False,
                    "qualifiers": {
                        "label": child_row["name"],
                        "database": child_row["db"],
                        "identity": f"{child_row['identity']:.1f}",
                        "match_length": f"{child_row['match']:.1f}",
                        "fragment": str(is_fragment),
                        "note": f"{child_row['db']}:{child_row['sseqid']}",
                    },
                }
            )

        records.append(
            {
                "key": f"{database}:{sequence_id}",
                "parent": {
                    "db": str(database),
                    "sseqid": str(sequence_id),
                    "name": str(name),
                    "type": str(feature_type),
                    "length": length,
                },
                "children": child_rows,
                "parsed": {
                    "name": str(sequence_id),
                    "length": length,
                    "topology": "linear",
                    "features": features,
                },
            }
        )
    return records


INDEX_HTML = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Nested feature figures — linear</title>
  <link rel="stylesheet" href="vendor/css/genes-d3.css">
  <link rel="stylesheet" href="viewer.css">
</head>
<body>
  <header class="masthead">
    <div>
      <p class="eyebrow">pLannotate database audit</p>
      <h1>Nested feature figures</h1>
      <p>Nested annotations · linear sequences · rendered with genes-d3</p>
    </div>
    <div class="summary" aria-label="audit totals">
      <span><strong id="parentCount">—</strong> parent features</span>
      <span><strong id="hitCount">—</strong> nested hits</span>
    </div>
  </header>

  <nav class="status-tabs" aria-label="curation bins">
    <button type="button" data-status="all">All <span id="allCount">—</span></button>
    <button type="button" data-status="good">Good <span id="goodCount">—</span></button>
    <button type="button" data-status="bad">Bad <span id="badCount">—</span></button>
    <button type="button" data-status="review">Review <span id="reviewCount">—</span></button>
  </nav>

  <main class="shell">
    <aside class="browser-panel">
      <label for="search">Find a parent or nested feature</label>
      <input id="search" type="search" placeholder="T5 promoter, lac operator…">
      <p id="matchCount" class="match-count"></p>
      <nav id="featureList" aria-label="problematic source features"></nav>
    </aside>

    <section class="figure-panel">
      <div class="figure-head">
        <div>
          <p id="parentKey" class="eyebrow"></p>
          <h2 id="parentName"></h2>
          <p id="parentMeta"></p>
        </div>
        <div class="actions">
          <button id="previous" type="button">Previous</button>
          <button id="next" type="button">Next</button>
          <button id="download" type="button">Download SVG</button>
        </div>
      </div>

      <div id="ruler"></div>
      <div class="map-layout">
        <div class="plate">
          <div id="plotArea"><div id="plot"></div></div>
          <div id="focus"></div>
        </div>
        <aside class="legend-panel">
          <label class="fragment-toggle">
            <input id="fragments" type="checkbox" checked>
            Show fragment calls
          </label>
          <div id="legend"></div>
        </aside>
      </div>

      <section class="calls">
        <div class="calls-heading">
          <h3>Nested annotation calls</h3>
          <p id="binExplanation"></p>
        </div>
        <div class="table-wrap">
          <table>
            <thead><tr><th>Feature</th><th>Source</th><th>Type</th><th>Coordinates</th><th>Identity</th><th>Coverage</th><th>E-value</th><th>Evidence</th><th>Bin</th><th>Rule</th><th>Action</th><th>Rationale</th></tr></thead>
            <tbody id="callRows"></tbody>
          </table>
        </div>
      </section>
    </section>
  </main>

  <script src="https://d3js.org/d3.v7.min.js"></script>
  <script src="vendor/js/parser.js"></script>
  <script src="vendor/js/data.js"></script>
  <script src="vendor/js/plot-core.js"></script>
  <script src="vendor/js/plot-linear.js"></script>
  <script src="vendor/js/plot-circular.js"></script>
  <script src="vendor/js/plot-stylized.js"></script>
  <script src="vendor/js/plot.js"></script>
  <script src="vendor/js/export.js"></script>
  <script src="audit-data.js"></script>
  <script src="viewer.js"></script>
</body>
</html>
"""

VIEWER_CSS = """* { box-sizing: border-box; }
body { margin: 0; background: var(--paper); color: var(--ink); font-family: var(--body); }
button, input { font: inherit; }
.masthead { display: flex; justify-content: space-between; gap: 24px; align-items: end; padding: 28px 32px 22px; background: var(--ink); color: var(--paper-hi); }
.masthead h1 { margin: 2px 0 6px; font-family: var(--label); font-size: clamp(1.7rem, 3vw, 2.6rem); }
.masthead p { margin: 0; color: rgba(248,244,235,.72); }
.eyebrow { margin: 0; color: var(--bronze) !important; font-family: var(--mono); font-size: .76rem; letter-spacing: .08em; text-transform: uppercase; }
.summary { display: flex; gap: 10px; flex-wrap: wrap; justify-content: end; }
.summary span { display: grid; padding: 10px 14px; border-radius: var(--r-md); background: rgba(248,244,235,.09); font-size: .78rem; }
.summary strong { color: var(--paper-hi); font-size: 1.25rem; }
.status-tabs { display: flex; gap: 8px; padding: 14px 18px 0; max-width: 1500px; margin: 0 auto; }
.status-tabs button { display: flex; gap: 7px; align-items: center; background: var(--paper-hi); color: var(--ink); box-shadow: inset 0 0 0 1px rgba(111,123,103,.22); }
.status-tabs button:hover { background: var(--paper-mid); }
.status-tabs button.active { background: var(--ink); color: var(--paper-hi); box-shadow: none; }
.status-tabs span { min-width: 1.7em; padding: 1px 6px; border-radius: var(--r-pill); background: rgba(111,123,103,.14); text-align: center; font-family: var(--mono); font-size: .72rem; }
.status-tabs button.active span { background: rgba(248,244,235,.15); }
.shell { display: grid; grid-template-columns: minmax(240px, 300px) minmax(0, 1fr); gap: 18px; padding: 18px; max-width: 1500px; margin: 0 auto; }
.browser-panel, .figure-panel { min-width: 0; }
.browser-panel { position: sticky; top: 18px; height: calc(100vh - 36px); padding: 18px; border-radius: var(--r-lg); background: var(--paper-hi); }
.browser-panel label { display: block; margin-bottom: 7px; font-weight: 700; }
#search { width: 100%; border: 0; border-radius: var(--r-pill); padding: 10px 13px; background: var(--paper-mid); color: var(--ink); outline: 2px solid transparent; }
#search:focus { outline-color: var(--bronze); }
.match-count { margin: 9px 4px; color: var(--ink-soft); font-size: .78rem; }
#featureList { height: calc(100% - 88px); overflow: auto; display: grid; align-content: start; gap: 5px; padding-right: 4px; }
.feature-link { width: 100%; border: 0; border-radius: 13px; padding: 9px 11px; text-align: left; background: transparent; color: var(--ink); cursor: pointer; }
.feature-link:hover { background: var(--paper-mid); }
.feature-link.active { background: var(--sage); color: var(--paper-hi); }
.feature-link strong, .feature-link small { display: block; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.feature-link small { margin-top: 3px; opacity: .72; }
.figure-panel { display: grid; grid-template-columns: minmax(0, 1fr); gap: 14px; }
.figure-head, #ruler, .plate, .calls { background: var(--paper-hi); border-radius: var(--r-lg); }
.figure-head { display: flex; justify-content: space-between; gap: 20px; align-items: center; padding: 20px 24px; }
.figure-head h2 { margin: 3px 0; font-size: 1.45rem; }
.figure-head p { margin: 0; color: var(--ink-soft); }
.actions { display: flex; gap: 7px; flex-wrap: wrap; justify-content: end; }
button { border: 0; border-radius: var(--r-pill); padding: 8px 12px; background: var(--sage); color: var(--paper-hi); cursor: pointer; }
button:hover { background: var(--sage-hi); }
button:disabled { opacity: .4; cursor: default; }
#ruler { padding: 12px 18px; overflow: hidden; }
.map-layout { display: grid; grid-template-columns: minmax(0, 1fr) 210px; gap: 14px; }
.plate { padding: 22px; min-height: 230px; }
#plot { max-width: 900px; margin: 0 auto; }
.legend-panel { padding: 16px; border-radius: var(--r-lg); background: var(--ink); color: var(--paper-hi); }
.fragment-toggle { display: flex; align-items: center; gap: 8px; margin-bottom: 12px; font-size: .82rem; }
.calls { padding: 20px 24px; }
.calls-heading { display: flex; align-items: baseline; justify-content: space-between; gap: 18px; margin-bottom: 12px; }
.calls h3, .calls-heading p { margin: 0; }
.calls-heading p { color: var(--ink-soft); font-size: .8rem; text-align: right; }
.table-wrap { overflow-x: auto; }
table { width: 100%; border-collapse: collapse; font-size: .84rem; }
th, td { padding: 8px 10px; border-top: 1px solid rgba(111,123,103,.18); text-align: left; white-space: nowrap; }
th { color: var(--ink-soft); font-size: .72rem; text-transform: uppercase; letter-spacing: .05em; }
td:first-child { font-weight: 700; }
.fragment-pill { color: var(--danger); }
.status-pill { display: inline-block; padding: 3px 8px; border-radius: var(--r-pill); font-weight: 700; text-transform: capitalize; }
.status-pill.good { background: rgba(80,125,80,.16); color: #356239; }
.status-pill.bad { background: rgba(171,70,53,.14); color: var(--danger); }
.status-pill.review { background: rgba(179,129,55,.18); color: #7b541d; }
.rationale { max-width: 38rem; white-space: normal; line-height: 1.35; }
@media (max-width: 920px) { .shell { grid-template-columns: 1fr; } .browser-panel { position: static; height: 300px; } .map-layout { grid-template-columns: 1fr; } }
@media (max-width: 620px) { .masthead, .figure-head, .calls-heading { align-items: start; flex-direction: column; } .status-tabs { overflow-x: auto; } .shell { padding: 10px; } .actions { justify-content: start; } .calls-heading p { text-align: left; } }
"""

VIEWER_JS = """(() => {
  const all = window.AUDIT_FIGURES || [];
  const validStatuses = new Set(["all", "good", "bad", "review"]);
  const requestedStatus = new URLSearchParams(location.search).get("status");
  let activeStatus = validStatuses.has(requestedStatus) ? requestedStatus : "all";
  let visible = [];
  let active = null;

  const byId = (id) => document.getElementById(id);
  const list = byId("featureList");
  const search = byId("search");
  const fragments = byId("fragments");
  const explanations = {
    all: "All rule outcomes",
    good: "Keep genuine nested biology or provenance; follow any record-correction action shown",
    bad: "Suppress false or unsupported calls, or repair a broken source parent",
    review: "No unresolved calls under the current formal policy",
  };

  function childrenFor(record) {
    return activeStatus === "all"
      ? record.children
      : record.children.filter((child) => child.status === activeStatus);
  }

  function searchable(record) {
    return [record.key, record.parent.name, record.parent.type]
      .concat(childrenFor(record).flatMap((child) => [child.name, child.sseqid, child.type, child.db, child.action, child.rationale]))
      .join(" ").toLowerCase();
  }

  function escapeHTML(value) {
    const span = document.createElement("span");
    span.textContent = String(value);
    return span.innerHTML;
  }

  function updateLocation() {
    const url = new URL(location.href);
    if (activeStatus === "all") url.searchParams.delete("status");
    else url.searchParams.set("status", activeStatus);
    url.hash = active ? encodeURIComponent(active.key) : "";
    history.replaceState(null, "", url);
  }

  function rebuildList() {
    list.replaceChildren();
    visible.forEach((record) => {
      const count = childrenFor(record).length;
      const button = document.createElement("button");
      button.type = "button";
      button.className = "feature-link" + (record.key === active?.key ? " active" : "");
      button.innerHTML = `<strong>${escapeHTML(record.parent.name)}</strong><small>${escapeHTML(record.key)} · ${count} hit${count === 1 ? "" : "s"}</small>`;
      button.addEventListener("click", () => select(record));
      list.appendChild(button);
    });
    byId("matchCount").textContent = `${visible.length.toLocaleString()} matching parent features`;
  }

  function parsedFor(record) {
    const selected = new Set(childrenFor(record));
    return {
      ...record.parsed,
      features: [record.parsed.features[0]].concat(
        record.parsed.features.slice(1).filter((feature, index) => selected.has(record.children[index]))
      ),
    };
  }

  function select(record, updateURL = true) {
    if (!record) return;
    active = record;
    const children = childrenFor(record);
    const parsed = parsedFor(record);
    byId("parentKey").textContent = record.key;
    byId("parentName").textContent = record.parent.name;
    byId("parentMeta").textContent = `${record.parent.type} · ${record.parent.length.toLocaleString()} bp · ${children.length} ${activeStatus === "all" ? "nested" : activeStatus} call${children.length === 1 ? "" : "s"}`;

    const result = PLOT.render(byId("plot"), parsed, {
      linear: true,
      stylized: false,
      fragments: fragments.checked,
    });
    PLOT_CORE.buildRuler(result.feats, parsed.length);
    PLOT_CORE.buildLegend(result.feats);
    renderTable(children);
    rebuildList();

    const index = visible.findIndex((item) => item.key === record.key);
    byId("previous").disabled = index <= 0;
    byId("next").disabled = index < 0 || index >= visible.length - 1;
    byId("download").disabled = false;
    if (updateURL) updateLocation();
    document.body.dataset.ready = "true";
  }

  function renderTable(children) {
    const body = byId("callRows");
    body.replaceChildren();
    children.forEach((child) => {
      const row = document.createElement("tr");
      const evalue = Number.isFinite(child.evalue) ? child.evalue.toExponential(1) : "—";
      row.innerHTML = `<td>${escapeHTML(child.name)}</td><td><code>${escapeHTML(child.db + ":" + child.sseqid)}</code></td><td>${escapeHTML(child.type)}</td><td>${child.start.toLocaleString()}–${child.end.toLocaleString()}</td><td>${child.identity.toFixed(1)}%</td><td>${child.match.toFixed(1)}%</td><td><code>${evalue}</code></td><td class="${child.fragment ? "fragment-pill" : ""}">${child.fragment ? "fragment" : "whole"}</td><td><span class="status-pill ${escapeHTML(child.status)}">${escapeHTML(child.status)}</span></td><td><code>${escapeHTML(child.decision_source)}</code></td><td><code>${escapeHTML(child.action)}</code></td><td class="rationale">${escapeHTML(child.rationale)}</td>`;
      body.appendChild(row);
    });
  }

  function preferredRecord() {
    const requested = decodeURIComponent(location.hash.slice(1));
    const preferredKeys = {
      all: "snapgene:MMLV_Psi",
      good: "snapgene:MMLV_Psi",
      bad: "snapgene:OpIE_2_promoter",
      review: "snapgene:CMV_intron_(3)",
    };
    const preferredKey = preferredKeys[activeStatus];
    return visible.find((record) => record.key === requested)
      || visible.find((record) => record.key === active?.key)
      || visible.find((record) => record.key === preferredKey)
      || visible[0];
  }

  function applyFilters() {
    const query = search.value.trim().toLowerCase();
    visible = all.filter((record) => childrenFor(record).length > 0)
      .filter((record) => !query || searchable(record).includes(query));
    const chosen = preferredRecord();
    active = chosen || null;
    byId("parentCount").textContent = visible.length.toLocaleString();
    byId("hitCount").textContent = visible.reduce((total, record) => total + childrenFor(record).length, 0).toLocaleString();
    rebuildList();
    if (chosen) {
      select(chosen);
    } else {
      byId("parentKey").textContent = "NO MATCHING CALLS";
      byId("parentName").textContent = activeStatus === "review" ? "Nothing requires review" : "No matching features";
      byId("parentMeta").textContent = activeStatus === "review"
        ? "Every current call is resolved by the formal policy."
        : "Clear the search or choose another curation bin.";
      ["plot", "ruler", "focus", "legend", "callRows"].forEach((id) => byId(id).replaceChildren());
      byId("previous").disabled = true;
      byId("next").disabled = true;
      byId("download").disabled = true;
      document.body.dataset.ready = "true";
    }
  }

  function setStatus(status) {
    activeStatus = status;
    document.querySelectorAll("[data-status]").forEach((button) => {
      button.classList.toggle("active", button.dataset.status === activeStatus);
      button.setAttribute("aria-pressed", String(button.dataset.status === activeStatus));
    });
    byId("binExplanation").textContent = explanations[activeStatus];
    applyFilters();
  }

  function step(delta) {
    const index = visible.findIndex((item) => item.key === active?.key);
    select(visible[index + delta]);
  }

  const counts = all.flatMap((record) => record.children).reduce((result, child) => {
    result[child.status] += 1;
    return result;
  }, { good: 0, bad: 0, review: 0 });
  byId("allCount").textContent = (counts.good + counts.bad + counts.review).toLocaleString();
  byId("goodCount").textContent = counts.good.toLocaleString();
  byId("badCount").textContent = counts.bad.toLocaleString();
  byId("reviewCount").textContent = counts.review.toLocaleString();

  document.querySelectorAll("[data-status]").forEach((button) => {
    button.addEventListener("click", () => setStatus(button.dataset.status));
  });
  search.addEventListener("input", applyFilters);
  fragments.addEventListener("change", () => PLOT.setFragments(fragments.checked));
  byId("previous").addEventListener("click", () => step(-1));
  byId("next").addEventListener("click", () => step(1));
  byId("download").addEventListener("click", () => {
    const svg = byId("plot").querySelector("svg");
    if (!svg || !active) return;
    const exported = EXPORT.svgText(svg);
    EXPORT.save(new Blob([exported.text], { type: "image/svg+xml" }), `${EXPORT.slug(active.key)}-linear.svg`);
  });

  setStatus(activeStatus);
})();
"""

README = """# Nested feature figures

Open `index.html` to browse all source features with nested annotations. Every map
uses genes-d3's linear renderer and the pLannotate nested-feature audit results.

The viewer applies the formal curation policy and exposes four tabs: all, good, bad,
and review. Good calls stay. Bad calls are either suppression candidates or broken
source-parent records marked for replacement. Whole children stay by default, so the
review bin is empty for the current audit. Every row shows its suggested action,
E-value, rule source, and rationale.

The `good/`, `bad/`, and `review/` directories open directly into the corresponding
filter. The good bin opens on MMLV Ψ, the bad bin on the wrong-frame OpIE2/AC152 call,
and review shows an explicit resolved-empty state. It supports searching, showing or
hiding fragment calls, and exporting the active linear map as SVG.

Regenerate from the repository root with:

```bash
python tools/nested_feature_figures.py --genes-d3 ~/projects/genes-d3
```

The `vendor/` files are copied from the explicitly supplied local genes-d3 checkout.
This viewer is a local inspection artifact and is excluded from the pLannotate source
distribution.
"""

BIN_INDEX = """<!doctype html>
<html lang="en">
<head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>{title} nested calls</title></head>
<body><p>Opening the {status} curation bin…</p><script>
const target = new URL("../index.html", location.href);
target.searchParams.set("status", "{status}");
target.hash = location.hash;
location.replace(target);
</script></body>
</html>
"""


def build_viewer(
    input_path: Path,
    output: Path,
    genes_d3: Path,
    decisions_path: Path | None = None,
) -> int:
    """Generate data, viewer files, and a pinned copy of genes-d3's assets."""
    report = pd.read_csv(input_path)
    required_columns = {
        "parent_db",
        "parent_sseqid",
        "parent_name",
        "parent_type",
        "parent_length",
        "nested_db",
        "nested_sseqid",
        "nested_name",
        "nested_type",
        "nested_strand",
        "nested_start",
        "nested_end",
        "percent_identity",
        "percent_match",
        "evalue",
        "fragment",
    }
    missing = required_columns - set(report.columns)
    if missing:
        raise ValueError(f"Audit CSV is missing columns: {', '.join(sorted(missing))}")

    report = classify_report(report)
    records = build_records(report)
    output.mkdir(parents=True, exist_ok=True)
    (output / "index.html").write_text(INDEX_HTML)
    (output / "viewer.css").write_text(VIEWER_CSS)
    (output / "viewer.js").write_text(VIEWER_JS)
    (output / "README.md").write_text(README)
    data = json.dumps(records, ensure_ascii=False, separators=(",", ":"))
    (output / "audit-data.js").write_text(f"window.AUDIT_FIGURES={data};\n")
    report.to_csv(output / "decisions.csv", index=False)
    if decisions_path is not None:
        decisions_path.parent.mkdir(parents=True, exist_ok=True)
        report.to_csv(decisions_path, index=False)

    for status in ("good", "bad", "review"):
        bin_directory = output / status
        bin_directory.mkdir(exist_ok=True)
        (bin_directory / "index.html").write_text(
            BIN_INDEX.format(title=status.title(), status=status)
        )

    vendor = output / "vendor"
    for relative_path in GENES_D3_ASSETS:
        source = genes_d3 / relative_path
        if not source.is_file():
            raise FileNotFoundError(f"Missing genes-d3 asset: {source}")
        destination = vendor / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    return len(records)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--genes-d3",
        type=Path,
        required=True,
        help="path to the local genes-d3 checkout whose renderer assets are copied",
    )
    parser.add_argument("--decisions", type=Path, default=DEFAULT_DECISIONS)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    count = build_viewer(args.input, args.output, args.genes_d3, args.decisions)
    print(f"Generated {count} linear nested-feature figures in {args.output}")


if __name__ == "__main__":
    main()
